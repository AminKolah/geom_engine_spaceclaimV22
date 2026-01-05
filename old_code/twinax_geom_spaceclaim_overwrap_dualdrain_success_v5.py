# Python Script, API Version = V22
# Python Script for SpaceClaim V22
import math
from SpaceClaim.Api.V22 import *
from SpaceClaim.Api.V22.Geometry import *
from SpaceClaim.Api.V22.Modeler import *

ClearAll()

# -----------------------------
# 1. PARAMETERS / INPUTS
# -----------------------------
def get_param(name, default_val):
    try: return float(Parameters[name])
    except: return default_val

# Logic options
filler_opt   = int(get_param('filler_option', 1))
drain_opt    = int(get_param('has_drains', 0))
is_a_doublet = int(get_param('is_a_doublet', 1)) 
is_elliptic  = int(get_param('is_elliptic', 1))
h_mix        = get_param('mixing_factor', 0.7)

# filler_mode logic
filler_mode = "shell" if filler_opt == 1 else "fill"

# Conductor and Core
D_cond    = get_param('diam_cond', 1.0)
C2C       = get_param('c2c', 3.0)
D_core    = get_param('diam_core', 2.8)

# Shield (Double-D)
t_shield   = get_param('t_shield', 0.1)
W_outer    = get_param('W_outer', 8.0)
H_outer    = get_param('H_outer', 5.0)

# Drains and Overwrap
D_drain    = get_param('diam_drain', 1.0)
t_overwrap = get_param('t_overwrap', 0.1)
L_extrude  = get_param('length_extrude', 50.0)

# Derived values
r_cond, r_core, r_drain = D_cond/2.0, D_core/2.0, D_drain/2.0
dx_cond = C2C / 2.0

# Shield inner dimensions
W_in = W_outer - 2.0 * t_shield
H_in = H_outer - 2.0 * t_shield
R_in = H_in / 2.0
dx_in = (W_in - H_in) / 2.0
R_shield = H_outer / 2.0
dx_shield = (W_outer - H_outer) / 2.0

# -----------------------------
# 2. HELPER FUNCTIONS
# -----------------------------
def set_sketch_plane_xy():
    ViewHelper.SetSketchPlane(Plane.PlaneXY)

def solidify_sketch():
    ViewHelper.SetViewMode(InteractionMode.Solid)

def sketch_circle(cx, cy, r):
    return SketchCircle.Create(Point2D.Create(MM(cx), MM(cy)), MM(r)).CreatedCurves

def sketch_doubleD(dx, R):
    c1 = SketchArc.CreateSweepArc(Point2D.Create(MM(dx), 0), Point2D.Create(MM(dx), MM(R)), Point2D.Create(MM(dx), MM(-R)), True).CreatedCurves
    c2 = SketchArc.CreateSweepArc(Point2D.Create(MM(-dx), 0), Point2D.Create(MM(-dx), MM(R)), Point2D.Create(MM(-dx), MM(-R)), False).CreatedCurves
    l1 = SketchLine.Create(Point2D.Create(MM(-dx), MM(R)), Point2D.Create(MM(dx), MM(R))).CreatedCurves
    l2 = SketchLine.Create(Point2D.Create(MM(-dx), MM(-R)), Point2D.Create(MM(dx), MM(-R))).CreatedCurves
    return list(c1) + list(c2) + list(l1) + list(l2)

def P2(x, y):
    return Point2D.Create(MM(x), MM(y))

def P2(x, y):
    return Point2D.Create(MM(x), MM(y))

def sketch_elliptic(W_val, H_val, h, n_ellipse=80):
    """
    Draws the mixed profile WITHOUT full side circles:
      - ellipse-only segment for |x| <= xc (top+bottom, polyline approx)
      - right/left circular end-cap ARCS only (no full circles)
    """

    D = W_val / 2.0
    H = H_val / 2.0

    # ---- Definitions from your slide ----
    C0 = D - H
    xc = C0 + h * H

    if not (0.0 < xc < D):
        raise Exception("Need 0 < xc < D. Got xc=%.6f, D=%.6f" % (xc, D))

    denom = (2.0 * D - xc)
    if abs(denom) < 1e-12:
        raise Exception("Invalid: 2D - xc ≈ 0 causes division by zero in rs.")

    rs = (H**2 + (D - xc) * D) / denom
    if (xc - D + rs) <= 0:
        raise Exception("Invalid: xc - D + rs must be > 0 for rx.")

    rx = math.sqrt((xc * H**2) / (xc - D + rs))
    if rs <= 0 or rx <= 0:
        raise Exception("Invalid geometry: rs or rx <= 0 (rs=%.6f, rx=%.6f)" % (rs, rx))

    # y at x=xc on ellipse
    arg = 1.0 - (xc * xc) / (rx * rx)
    if arg <= 0:
        raise Exception("Ellipse does not reach x=xc (arg=%.6f). Check h/W/H." % arg)
    yc = H * math.sqrt(arg)

    # ---- 1) Ellipse portion (|x|<=xc), approximated by line segments ----
    xs = [(-xc + (2.0 * xc) * i / float(n_ellipse)) for i in range(n_ellipse + 1)]

    # top
    pts_top = [P2(x,  H * math.sqrt(max(0.0, 1.0 - (x*x)/(rx*rx)))) for x in xs]
    for i in range(len(pts_top) - 1):
        SketchLine.Create(pts_top[i], pts_top[i+1])

    # bottom (reverse)
    pts_bot = [P2(x, -H * math.sqrt(max(0.0, 1.0 - (x*x)/(rx*rx)))) for x in reversed(xs)]
    for i in range(len(pts_bot) - 1):
        SketchLine.Create(pts_bot[i], pts_bot[i+1])

    # ---- 2) End-cap arcs ONLY (no full circles) ----
    # Right circle center at (+a, 0), where a = D - rs
    a = D - rs

    # angle of join point on right circle (computed from the join point location)
    # vector from center to join is (xc - a, yc)
    phi = math.atan2(yc, (xc - a))

    # Build points exactly on circle (robust vs tiny mismatch)
    # Right side join points:
    xr = a + rs * math.cos(phi)
    yr =     rs * math.sin(phi)

    start_r = P2(xr, +yr)
    end_r   = P2(xr, -yr)

    # Right cap arc: from top -> bottom going through (D,0) (rightmost), i.e. clockwise
    SketchArc.CreateSweepArc(P2(+a, 0.0), start_r, end_r, True)

    # Left circle center at (-a, 0)
    # Join points mirrored:
    xl = -a - rs * math.cos(phi)
    yl =      rs * math.sin(phi)

    start_l = P2(xl, +yl)
    end_l   = P2(xl, -yl)

    # Left cap arc: top -> bottom going through (-D,0) (leftmost), i.e. counterclockwise
    SketchArc.CreateSweepArc(P2(-a, 0.0), start_l, end_l, False)

    print("Elliptic mix (arcs-only) sketched:",
          "xc=%.4f rx=%.4f rs=%.4f yc=%.4f" % (xc, rx, rs, yc))

def sketch_profile(dx, R, W_val, H_val):
    if is_elliptic:
        sketch_elliptic(W_val, H_val, h_mix)
    else:
        sketch_doubleD(dx, R)

def extrude_and_name(name, length_mm):
    solidify_sketch()
    temp_body = GetRootPart().Bodies[-1]
    
    # Instead of picking [0] or the largest, select ALL faces
    # This treats the overlapping ellipse and circles as ONE object
    all_faces = FaceSelection.Create(temp_body.Faces)
    
    options = ExtrudeFaceOptions()
    options.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(all_faces, Direction.DirZ, MM(length_mm), options)
    
    new_body = res.CreatedBodies[0]
    new_body.SetName(name)
    # This step removes internal "imprint" lines to make it neat
    # Part.Check(new_body) or similar can be used here.

# -----------------------------
# 3. BUILD GEOMETRY
# -----------------------------

# (1) Conductors
set_sketch_plane_xy()
sketch_circle(-dx_cond, 0, r_cond)
extrude_and_name("conductor[1]", L_extrude)

set_sketch_plane_xy()
sketch_circle(dx_cond, 0, r_cond)
extrude_and_name("conductor[2]", L_extrude)

# (2) Cores (Insulation)
if not (is_a_doublet):
    set_sketch_plane_xy()
    sketch_circle(-dx_cond, 0, r_core)
    sketch_circle(-dx_cond, 0, r_cond)
    extrude_and_name("single_core[1]", L_extrude, True)

    set_sketch_plane_xy()
    sketch_circle(dx_cond, 0, r_core)
    sketch_circle(dx_cond, 0, r_cond)
    extrude_and_name("single_core[2]", L_extrude, True)

# (3) Filler / Doublet
set_sketch_plane_xy()
if is_a_doublet and filler_opt == 1:
    sketch_profile(dx_in, R_in, W_in, H_in)
    solidify_sketch()
    sketch_circle(-dx_cond, 0, r_cond)
    sketch_circle(dx_cond, 0, r_cond)
    extrude_and_name("Second_Extrusion", L_extrude)
