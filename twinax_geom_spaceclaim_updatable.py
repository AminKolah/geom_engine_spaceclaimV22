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
filler_opt   = int(get_param('filler_option', 0))
drain_opt    = int(get_param('has_drains', 1))
is_a_doublet = int(get_param('is_a_doublet', 0)) 
is_elliptic  = int(get_param('is_elliptic', 0))
h_mix        = get_param('mixing_factor', 0.9)

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

    # ---- 1) Ellipse portion (|x|<=xc) using smooth Splines ----
    xs = [(-xc + (2.0 * xc) * i / float(n_ellipse)) for i in range(n_ellipse + 1)]

    # TOP Spline
    # Collect points into a List[Point2D]
    pts_top = [P2(x, H * math.sqrt(max(0.0, 1.0 - (x*x)/(rx*rx)))) for x in xs]
    # Create a single smooth curve through all points
    SketchNurbs.CreateFrom2DPoints(False, pts_top) 

    # BOTTOM Spline
    pts_bot = [P2(x, -H * math.sqrt(max(0.0, 1.0 - (x*x)/(rx*rx)))) for x in reversed(xs)]
    SketchNurbs.CreateFrom2DPoints(False, pts_bot)

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

def sketch_wrap_loop_with_drains(offset):
    if not is_elliptic:
        d = abs(x_drain - dx_shield) 
        R, r = R_shield, r_drain
        theta = math.asin((R - r) / d)
        Ro, ro = R + offset, r + offset
        tx_S, ty_S = Ro * math.sin(theta), Ro * math.cos(theta)
        tx_D, ty_D = ro * math.sin(theta), ro * math.cos(theta)
        SketchArc.CreateSweepArc(P2(x_drain, 0), P2(x_drain + tx_D, ty_D), P2(x_drain + tx_D, -ty_D), True)
        SketchArc.CreateSweepArc(P2(-x_drain, 0), P2(-x_drain - tx_D, ty_D), P2(-x_drain - tx_D, -ty_D), False)
        SketchLine.Create(P2(x_drain + tx_D, ty_D), P2(dx_shield + tx_S, ty_S))
        SketchLine.Create(P2(x_drain + tx_D, -ty_D), P2(dx_shield + tx_S, -ty_S))
        SketchLine.Create(P2(-x_drain - tx_D, ty_D), P2(-dx_shield - tx_S, ty_S))
        SketchLine.Create(P2(-x_drain - tx_D, -ty_D), P2(-dx_shield - tx_S, -ty_S))
        SketchArc.CreateSweepArc(P2(dx_shield, 0), P2(dx_shield + tx_S, ty_S), P2(dx_shield, Ro), False)
        SketchArc.CreateSweepArc(P2(dx_shield, 0), P2(dx_shield + tx_S, -ty_S), P2(dx_shield, -Ro), True)
        SketchArc.CreateSweepArc(P2(-dx_shield, 0), P2(-dx_shield, Ro), P2(-dx_shield - tx_S, ty_S), False)
        SketchArc.CreateSweepArc(P2(-dx_shield, 0), P2(-dx_shield, -Ro), P2(-dx_shield - tx_S, -ty_S), True)
        SketchLine.Create(P2(-dx_shield, Ro), P2(dx_shield, Ro))
        SketchLine.Create(P2(-dx_shield, -Ro), P2(dx_shield, -Ro))
    else:  
        ro = r_drain + offset
        # Parameters for the overwrap ellipse
        Ho = (H_outer / 2.0) + offset
        Do = (W_outer / 2.0) + offset
        
        # Calculate rx for the overwrap (using your mixing logic)
        C0_o = Do - Ho
        xc_o = C0_o + h_mix * Ho
        rs_o = (Ho**2 + (Do - xc_o) * Do) / (2.0 * Do - xc_o)
        rx_o = math.sqrt((xc_o * Ho**2) / (xc_o - Do + rs_o))

        # 1. Calculate the 'Shoulder' of the drain where the line will land
        # We want a line tangent to both the ellipse and the circle.
        # A good approximation for the 'breakaway' angle theta:
        d = abs(x_drain) # Distance to drain center
        # Theta is the angle of the tangent line
        theta = math.asin((Ho - ro) / d) 
        
        # 2. Tangent point on the Circle (Drain)
        tx_D = ro * math.sin(theta)
        ty_D = ro * math.cos(theta)
        
        # 3. Tangent point on the Ellipse
        # For a slope m = -tan(theta), the tangent x on an ellipse is:
        # x = (a^2 * m) / sqrt(a^2 * m^2 + b^2)
        m = math.tan(theta)
        xt = (rx_o**2 * m) / math.sqrt(rx_o**2 * m**2 + Ho**2)
        yt = Ho * math.sqrt(max(0, 1.0 - (xt**2 / rx_o**2)))

        # 4. DRAWING
        # Drain Arcs (Outer part)
        # Right Drain
        SketchArc.CreateSweepArc(P2(x_drain, 0), P2(x_drain + tx_D, ty_D), P2(x_drain + tx_D, -ty_D), True)
        # Left Drain
        SketchArc.CreateSweepArc(P2(-x_drain, 0), P2(-x_drain - tx_D, ty_D), P2(-x_drain - tx_D, -ty_D), False)
        
        # Tangent Bridge Lines
        # Right side
        SketchLine.Create(P2(x_drain + tx_D, ty_D), P2(xt, yt))
        SketchLine.Create(P2(x_drain + tx_D, -ty_D), P2(xt, -yt))
        # Left side
        SketchLine.Create(P2(-x_drain - tx_D, ty_D), P2(-xt, yt))
        SketchLine.Create(P2(-x_drain - tx_D, -ty_D), P2(-xt, -yt))
        
        # 5. Top/Bottom Elliptic Curves (between -xt and xt)
        n_step = 20
        xs = [(-xt + (2.0 * xt) * i / float(n_step)) for i in range(n_step + 1)]
        
        pts_top = [P2(x, Ho * math.sqrt(max(0.0, 1.0 - (x*x)/(rx_o*rx_o)))) for x in xs]
        pts_bot = [P2(x, -Ho * math.sqrt(max(0.0, 1.0 - (x*x)/(rx_o*rx_o)))) for x in reversed(xs)]
        
        SketchNurbs.CreateFrom2DPoints(False, pts_top)
        SketchNurbs.CreateFrom2DPoints(False, pts_bot)
        
def extrude_and_name(name, length_mm, pick_largest=False):
    solidify_sketch()
    temp_body = GetRootPart().Bodies[-1]
    if pick_largest:
        target_face = sorted(temp_body.Faces, key=lambda f: f.Area, reverse=True)[0]
    else:
        target_face = temp_body.Faces[0]
    options = ExtrudeFaceOptions()
    options.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(target_face), Direction.DirZ, MM(length_mm), options)
    new_body = res.CreatedBodies[0]
    new_body.SetName(name)
    return new_body

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
    sketch_circle(-dx_cond, 0, r_cond)
    sketch_circle(dx_cond, 0, r_cond)
    extrude_and_name("Second_Extrusion", L_extrude, True)

elif filler_mode == "fill":
    sketch_profile(dx_in, R_in, W_in, H_in)
    sketch_circle(-dx_cond, 0, r_core)
    sketch_circle(dx_cond, 0, r_core)
    extrude_and_name("Second_Extrusion", L_extrude, True)

elif filler_mode == "shell":
    t_filler_shell = (dx_in + R_in) - (dx_cond + r_core)
    sketch_profile(dx_in, R_in, W_in, H_in)
    
    W_fill_inner = W_in - 2.0 * t_filler_shell
    H_fill_inner = H_in - 2.0 * t_filler_shell
    R_fill_inner = R_in - t_filler_shell
    
    if H_fill_inner > 0:
        sketch_profile(dx_in, R_fill_inner, W_fill_inner, H_fill_inner)
    
    sketch_circle(-dx_cond, 0, r_core)
    sketch_circle(dx_cond, 0, r_core)
    extrude_and_name("Second_Extrusion", L_extrude, True)

# (4) Shield
set_sketch_plane_xy()
sketch_profile(dx_shield, R_shield, W_outer, H_outer)
extrude_and_name("Shield", L_extrude, True)

# (5) Drains
x_drain = dx_shield + R_shield + r_drain
if drain_opt:
    set_sketch_plane_xy()
    sketch_circle(-x_drain, 0, r_drain)
    extrude_and_name("drain[1]", L_extrude)
    set_sketch_plane_xy()
    sketch_circle(x_drain, 0, r_drain)
    extrude_and_name("drain[2]", L_extrude)

# (6) Overwrap
if drain_opt:
    # Existing drain-wrap logic relies on Double-D tangency; 
    # This remains as previously defined
    set_sketch_plane_xy()
    sketch_wrap_loop_with_drains(0)
    sketch_wrap_loop_with_drains(t_overwrap)
    extrude_and_name("Overwrap", L_extrude, True)
else:
    set_sketch_plane_xy()
    # Fixed argument mismatch here
    sketch_profile(dx_shield, R_shield + t_overwrap, W_outer + 2*t_overwrap, H_outer + 2*t_overwrap)
    extrude_and_name("Overwrap", L_extrude, True)

# -----------------------------
# 4. CLEANUP
# -----------------------------
#surfaces = DataModel.GetObjectsByName("Surface")
#if surfaces.Count > 0: Selection.Create(surfaces).Delete()
#Selection.Create(GetRootPart().GetAllCurves()).Delete()

print("Assembly Created. Doublet: " + str(bool(is_a_doublet)) + " Elliptic: " + str(bool(is_elliptic)))
