

# Python Script, API Version = V252
import math
from SpaceClaim.Api.V252 import *
from SpaceClaim.Api.V252.Geometry import *
from SpaceClaim.Api.V252.Modeler import *
import sys
ClearAll()

# -----------------------------
# 1. PARAMETERS / INPUTS
# -----------------------------

def get_param(name, default_val):
    try: return float(Parameters[name])
    except: return default_val


run_benchmark_rig = False  # Set to True to enable 3-point rig creation
run_stiffness = False  # set False to silence

# Logic options using the Parameters check pattern

filler_opt   = int(Parameters.filler_option) if hasattr(Parameters, 'filler_option') else 0
drain_opt    = float(Parameters.has_drains) if hasattr(Parameters, 'has_drains') else 0.0
is_a_doublet = int(Parameters.is_a_doublet) if hasattr(Parameters, 'is_a_doublet') else 0
is_elliptic  = int(Parameters.is_elliptic) if hasattr(Parameters, 'is_elliptic') else 1
h_mix        = float(Parameters.mixing_factor) if hasattr(Parameters, 'mixing_factor') else 0.7
bending_benchmark_opt = int(Parameters.bending_mode_option) if hasattr(Parameters, 'bending_mode_option') else 0



if bending_benchmark_opt == 0:
    bending_mode = "good_3point"
elif bending_benchmark_opt == 1:
    bending_mode = "bad_3point"
elif bending_benchmark_opt == 2:
    bending_mode = "good_2pulleys"
elif bending_benchmark_opt == 3:
    bending_mode = "bad_2pulleys"
else:
    raise Exception("Unknown bending_mode_option: %s" % bending_benchmark_opt)

# Conductor and Core

D_cond = float(Parameters.diam_cond) if hasattr(Parameters, 'diam_cond') else 1.0
core_overlap = float(Parameters.core_overlap) if hasattr(Parameters, 'core_overlap') else 0.01
if core_overlap < 0:
    raise Exception("core_overlap must be >= 0. Negative values would create a real gap between cores which we don't want.")
D_core = float(Parameters.diam_core) if hasattr(Parameters, 'diam_core') else 2.8
# representation choice only
merge_cores = int(Parameters.merge_cores) if hasattr(Parameters, 'merge_cores') else 0
core_mode = "merged" if merge_cores else "separate"


# filler_mode logic

filler_mode = "shell" if filler_opt == 1 else "fill"

min_overlap_fill = 0.001  # a min overlap
if filler_mode == "fill" and core_overlap < min_overlap_fill:
        raise Exception(
        "For filled second extrusion, cores need a positive overlap to avoid zero-thickness filler geometry. "
        "Got core_overlap=%.6f mm." % core_overlap)


# Shield (Double-D)

t_shield = float(Parameters.t_shield) if hasattr(Parameters, 't_shield') else 0.1
W_outer  = float(Parameters.W_outer) if hasattr(Parameters, 'W_outer') else 8.0
H_outer  = float(Parameters.H_outer) if hasattr(Parameters, 'H_outer') else 5.0

# Drains and Overwrap

D_drain    = float(Parameters.diam_drain) if hasattr(Parameters, 'diam_drain') else 1.0
t_overwrap = float(Parameters.t_overwrap) if hasattr(Parameters, 't_overwrap') else 0.1
L_extrude  = float(Parameters.length_extrude) if hasattr(Parameters, 'length_extrude') else 50.0

# Logic Options

filler_opt   = int(Parameters.filler_option) if hasattr(Parameters, 'filler_option') else 0
drain_opt    = int(Parameters.has_drains) if hasattr(Parameters, 'has_drains') else 0
is_a_doublet = int(Parameters.is_a_doublet) if hasattr(Parameters, 'is_a_doublet') else 0 
is_elliptic  = int(Parameters.is_elliptic) if hasattr(Parameters, 'is_elliptic') else 1
h_mix        = float(Parameters.mixing_factor) if hasattr(Parameters, 'mixing_factor') else 0.7

# Multilumen options

is_a_multilumen           = int(Parameters.is_a_multilumen) if hasattr(Parameters, 'is_a_multilumen') else 0
multilumen_shape_opt      = int(Parameters.multilumen_shape) if hasattr(Parameters, 'multilumen_shape') else 0
multilumen_wall_thickness = float(Parameters.multilumen_wall_thickness) if hasattr(Parameters, 'multilumen_wall_thickness') else 0.1
n_multilumen_cavity       = int(Parameters.n_multilumen_cavity) if hasattr(Parameters, 'n_multilumen_cavity') else 9

# multilumen shape mode logic

multilumen_shape_mode = "trapezoidal" if multilumen_shape_opt == 1 else "circular"

# Derived values

r_cond, r_core, r_drain = D_cond/2.0, D_core/2.0, D_drain/2.0
C2C = D_core - core_overlap
dx_cond = 0.5 * C2C

# Shield inner dimensions

W_in = W_outer - 2.0 * t_shield
H_in = H_outer - 2.0 * t_shield
R_in = H_in / 2.0
dx_in = (W_in - H_in) / 2.0
R_shield = H_outer / 2.0
dx_shield = (W_outer - H_outer) / 2.0


print("Driven C2C = %.6f mm" % C2C)
print("Driven dx_cond = %.6f mm" % dx_cond)



# -----------------------------
# 2. HELPER FUNCTIONS
# -----------------------------
def clear_all_sketch_curves():
    root = GetRootPart()
    try:
        crvs = list(root.Curves)
        if crvs:
            Delete.Execute(CurveSelection.Create(crvs))
    except:
        pass

def ClearAllSketchCurves():
    """
    Finds all curves in the active design that are not yet part of a 3D solid 
    and deletes them to provide a clean slate for the next sketch.
    """
    root = GetRootPart()
    # Collect all curves currently in the Part
    curves = list(root.Curves)
    if len(curves) > 0:
        # Create a selection and delete them
        Delete.Execute(CurveSelection.Create(curves))
        

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


    curves = [] # Create a list to hold all the curves


    # 1) TOP Spline

    pts_top = [P2(x, H * math.sqrt(max(0.0, 1.0 - (x*x)/(rx*rx)))) for x in xs]
    res_top = SketchNurbs.CreateFrom2DPoints(False, pts_top)
    curves.extend(res_top.CreatedCurves) # Add to list


    # BOTTOM Spline

    pts_bot = [P2(x, -H * math.sqrt(max(0.0, 1.0 - (x*x)/(rx*rx)))) for x in reversed(xs)]
    res_bot = SketchNurbs.CreateFrom2DPoints(False, pts_bot)
    curves.extend(res_bot.CreatedCurves) # Add to list


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

    res_r = SketchArc.CreateSweepArc(P2(+a, 0.0), start_r, end_r, True)
    curves.extend(res_r.CreatedCurves)

    # Left circle center at (-a, 0)
    # Join points mirrored:

    xl = -a - rs * math.cos(phi)
    yl =      rs * math.sin(phi)

    start_l = P2(xl, +yl)
    end_l   = P2(xl, -yl)

    # Left cap arc: top -> bottom going through (-D,0) (leftmost), i.e. counterclockwise

    res_l = SketchArc.CreateSweepArc(P2(-a, 0.0), start_l, end_l, False)
    curves.extend(res_l.CreatedCurves)


    print("Elliptic mix (arcs-only) sketched:",
          "xc=%.4f rx=%.4f rs=%.4f yc=%.4f" % (xc, rx, rs, yc))

    return curves # IMPORTANT: Return the list of curves

def sketch_trapezoid(cx, cy, angle_rad, h, w_inner, w_outer):
    """
    Draw a trapezoid lumen centered at (cx,cy), whose symmetry axis is along angle_rad.
    "Inner" side is toward -radial direction, "Outer" toward +radial direction.
      - height = h (radial direction)
      - w_inner = width at inner side
      - w_outer = width at outer side
    """
    
    ca = math.cos(angle_rad)
    sa = math.sin(angle_rad)

    # Local axes: radial u (outward), tangential v (CCW)
    ux, uy = ca, sa
    vx, vy = -sa, ca

    # Half-sizes
    hi = 0.5 * h
    wi = 0.5 * w_inner
    wo = 0.5 * w_outer

    # Four corners in (u,v) local coords:
    # inner face at u = -hi, outer face at u = +hi
    # inner width = w_inner, outer width = w_outer
    local = [
        (-hi, +wi),  # inner-top
        (+hi, +wo),  # outer-top
        (+hi, -wo),  # outer-bot
        (-hi, -wi),  # inner-bot
    ]

    pts = []
    for (u, v) in local:
        x = cx + u*ux + v*vx
        y = cy + u*uy + v*vy
        pts.append(P2(x, y))

    # Connect
    SketchLine.Create(pts[0], pts[1])
    SketchLine.Create(pts[1], pts[2])
    SketchLine.Create(pts[2], pts[3])
    SketchLine.Create(pts[3], pts[0])


def sketch_profile(dx, R, W_val, H_val):
    
    if is_elliptic:
        return sketch_elliptic(W_val, H_val, h_mix)
    else:
        return sketch_doubleD(dx, R)


def sketch_wrap_loop_with_drains(offset):
    curves = []

    if not is_elliptic:
        d = abs(x_drain - dx_shield)
        R, r = R_shield, r_drain
        theta = math.asin((R - r) / d)
        Ro, ro = R + offset, r + offset
        tx_S, ty_S = Ro * math.sin(theta), Ro * math.cos(theta)
        tx_D, ty_D = ro * math.sin(theta), ro * math.cos(theta)

        curves.extend(SketchArc.CreateSweepArc(P2(x_drain, 0), P2(x_drain + tx_D, ty_D), P2(x_drain + tx_D, -ty_D), True).CreatedCurves)
        curves.extend(SketchArc.CreateSweepArc(P2(-x_drain, 0), P2(-x_drain - tx_D, ty_D), P2(-x_drain - tx_D, -ty_D), False).CreatedCurves)

        curves.extend(SketchLine.Create(P2(x_drain + tx_D, ty_D), P2(dx_shield + tx_S, ty_S)).CreatedCurves)
        curves.extend(SketchLine.Create(P2(x_drain + tx_D, -ty_D), P2(dx_shield + tx_S, -ty_S)).CreatedCurves)
        curves.extend(SketchLine.Create(P2(-x_drain - tx_D, ty_D), P2(-dx_shield - tx_S, ty_S)).CreatedCurves)
        curves.extend(SketchLine.Create(P2(-x_drain - tx_D, -ty_D), P2(-dx_shield - tx_S, -ty_S)).CreatedCurves)

        curves.extend(SketchArc.CreateSweepArc(P2(dx_shield, 0), P2(dx_shield + tx_S, ty_S), P2(dx_shield, Ro), False).CreatedCurves)
        curves.extend(SketchArc.CreateSweepArc(P2(dx_shield, 0), P2(dx_shield + tx_S, -ty_S), P2(dx_shield, -Ro), True).CreatedCurves)
        curves.extend(SketchArc.CreateSweepArc(P2(-dx_shield, 0), P2(-dx_shield, Ro), P2(-dx_shield - tx_S, ty_S), False).CreatedCurves)
        curves.extend(SketchArc.CreateSweepArc(P2(-dx_shield, 0), P2(-dx_shield, -Ro), P2(-dx_shield - tx_S, -ty_S), True).CreatedCurves)

        curves.extend(SketchLine.Create(P2(-dx_shield, Ro), P2(dx_shield, Ro)).CreatedCurves)
        curves.extend(SketchLine.Create(P2(-dx_shield, -Ro), P2(dx_shield, -Ro)).CreatedCurves)

    else:
        curves = []
        ro = r_drain + offset
        Ho = (H_outer / 2.0) + offset
        Do = (W_outer / 2.0) + offset

        C0_o = Do - Ho
        xc_o = C0_o + h_mix * Ho
        rs_o = (Ho**2 + (Do - xc_o) * Do) / (2.0 * Do - xc_o)
        rx_o = math.sqrt((xc_o * Ho**2) / (xc_o - Do + rs_o))

        d = abs(x_drain)
        theta = math.asin((Ho - ro) / d)

        tx_D = ro * math.sin(theta)
        ty_D = ro * math.cos(theta)

        m = math.tan(theta)
        xt = (rx_o**2 * m) / math.sqrt(rx_o**2 * m**2 + Ho**2)
        yt = Ho * math.sqrt(max(0, 1.0 - (xt**2 / rx_o**2)))

        curves.extend(SketchArc.CreateSweepArc(P2(x_drain, 0), P2(x_drain + tx_D, ty_D), P2(x_drain + tx_D, -ty_D), True).CreatedCurves)
        curves.extend(SketchArc.CreateSweepArc(P2(-x_drain, 0), P2(-x_drain - tx_D, ty_D), P2(-x_drain - tx_D, -ty_D), False).CreatedCurves)

        curves.extend(SketchLine.Create(P2(x_drain + tx_D, ty_D), P2(xt, yt)).CreatedCurves)
        curves.extend(SketchLine.Create(P2(x_drain + tx_D, -ty_D), P2(xt, -yt)).CreatedCurves)
        curves.extend(SketchLine.Create(P2(-x_drain - tx_D, ty_D), P2(-xt, yt)).CreatedCurves)
        curves.extend(SketchLine.Create(P2(-x_drain - tx_D, -ty_D), P2(-xt, -yt)).CreatedCurves)

        n_step = 20
        xs = [(-xt + (2.0 * xt) * i / float(n_step)) for i in range(n_step + 1)]

        pts_top = [P2(x, Ho * math.sqrt(max(0.0, 1.0 - (x*x)/(rx_o*rx_o)))) for x in xs]
        pts_bot = [P2(x, -Ho * math.sqrt(max(0.0, 1.0 - (x*x)/(rx_o*rx_o)))) for x in reversed(xs)]

        curves.extend(SketchNurbs.CreateFrom2DPoints(False, pts_top).CreatedCurves)
        curves.extend(SketchNurbs.CreateFrom2DPoints(False, pts_bot).CreatedCurves)

    return curves



def sketch_wrap_loop_with_drains_old(offset):
    
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

def share_topology_pair(body_a, body_b, mode="share_only", verbose=True):
    """
    mode:
      - "share_only": sets SharedTopology only (HFSS-cleanest)
      - "share_and_imprint": sets SharedTopology + tries Imprint/Repair.Share (stronger, may fragment faces)
    """
    if body_a is None or body_b is None:
        return

    # --- 1) Set SharedTopology on BOTH parents (not just body_a) ---
    for b in [body_a, body_b]:
        try:
            comp = b.Parent
            if hasattr(comp, 'SharedTopology'):
                comp.SharedTopology = SharedTopology.Share
                if verbose:
                    print("Component Property set to 'Share' for:", b.Name)
        except Exception as e:
            if verbose:
                print("Failed to set SharedTopology property for", b.Name, ":", e)

    # --- 2) Optionally imprint ---
    if mode != "share_and_imprint":
        return

    try:
        selection = BodySelection.Create([body_a, body_b])
        options = ImprintOptions()
        Imprint.Execute(selection, options)
        if verbose:
            print("Physical Imprint successful between:", body_a.Name, "and", body_b.Name)
    except Exception as e_imprint:
        if verbose:
            print("Physical Imprint failed:", e_imprint)
        try:
            Share.Share(selection)
            if verbose:
                print("Repair.Share successful.")
        except:
            if verbose:
                print("All sharing methods failed.")

def _largest_xy_face(body):
    best = None
    bestA = -1.0
    for f in list(body.Faces):
        try:
            n = f.Plane.Normal
            # face normal near ±Z indicates XY face
            if abs(n.Z) > 0.95:
                if f.Area > bestA:
                    best = f
                    bestA = f.Area
        except:
            pass
    return best


def area_doubleD(W, H):
    R = H / 2.0
    dx = (W - H) / 2.0
    return 4.0 * dx * R + math.pi * R * R


def extrude_and_name(name, length_mm, pick_largest=True):
    root = GetRootPart()
    before = list(root.Bodies)

    solidify_sketch()

    after = list(root.Bodies)
    new_bodies = [b for b in after if b not in before]

    if not new_bodies:
        try:
            crvs = list(root.Curves)
            if not crvs:
                raise Exception("No sketch curves found to Fill for '%s'." % name)
            Fill.Execute(CurveSelection.Create(crvs))
        except Exception as e:
            raise Exception("Solidify produced NO sketch-fill body for '%s' and Fill failed: %s" % (name, e))

        after2 = list(root.Bodies)
        new_bodies = [b for b in after2 if b not in before]
        if not new_bodies:
            raise Exception("Solidify+Fill produced NO sketch-fill body for '%s'." % name)

    def best_xy_area(b):
        f = _largest_xy_face(b)
        if f is None:
            try:
                f = max(list(b.Faces), key=lambda ff: ff.Area)
            except:
                return -1.0
        try:
            return float(f.Area)
        except:
            return -1.0

    if pick_largest:
        temp_body = max(new_bodies, key=best_xy_area)
    else:
        bodies_sorted = sorted(new_bodies, key=best_xy_area)
        temp_body = bodies_sorted[len(bodies_sorted)//2]

    face = _largest_xy_face(temp_body)
    if face is None:
        face = max(list(temp_body.Faces), key=lambda f: f.Area)

    options = ExtrudeFaceOptions()
    options.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirZ, MM(length_mm), options)

    created = list(res.CreatedBodies)
    if not created:
        raise Exception("Extrude produced no result for '%s'." % name)

    new_body = created[0]
    new_body.SetName(name)

    try:
        Delete.Execute(BodySelection.Create(new_bodies))
    except:
        pass
    
    try:
        clear_all_sketch_curves()
    except:
        pass

    return new_body


def _get_named_selection_by_name(ns_name):
    try:
        root = GetRootPart()
        for ns in list(root.NamedSelections):
            try:
                if ns.GetName() == ns_name:
                    return ns
            except:
                if getattr(ns, "Name", "") == ns_name:
                    return ns
    except:
        pass
    return None

def _delete_named_selection_if_exists(ns_name):
    ns = _get_named_selection_by_name(ns_name)
    if ns is None:
        return
    try:
        Delete.Execute(Selection.Create(ns))
    except:
        try:
            ns.Delete()
        except:
            pass

def union_bodies(name, bodies):
    """
    Robust union/merge for SpaceClaim.
    - Tries Boolean.Unite first (most reliable for solids)
    - Falls back to Combine.Merge patterns used in some builds
    Returns the resulting body (best effort).
    """
    bodies = [b for b in bodies if b is not None]
    if len(bodies) < 2:
        return bodies[0] if bodies else None

    # Ensure they are in the root part's body list (or same component) if possible
    # (Not strictly required, but helps avoid some failures.)

    # --- Try Boolean.Unite (preferred) ---
    try:
        # Often expects a Selection, not BodySelection
        res = Boolean.Unite(Selection.Create(bodies))
        # Different return shapes across builds
        try:
            out = list(res.CreatedBodies)[0]
        except:
            try:
                out = res.Body
            except:
                out = bodies[0]
        try:
            out.SetName(name)
        except:
            pass
        return out
    except Exception as e1:
        print("Boolean.Unite failed:", e1)

    # --- Try Combine.Merge variants ---
    try:
        res = Combine.Merge(Selection.Create(bodies))
        # Sometimes you get MergedBody, sometimes CreatedBodies
        try:
            out = res.MergedBody
        except:
            try:
                out = list(res.CreatedBodies)[0]
            except:
                out = bodies[0]
        try:
            out.SetName(name)
        except:
            pass
        return out
    except Exception as e2:
        print("Combine.Merge failed:", e2)

    # If we got here, union failed
    raise Exception("Union failed. Bodies may not intersect/touch, or are invalid (sheet, separate components, etc.).")


def _face_bbox_xy(face):
    bb = face.Shape.GetBoundingBox(Matrix.Identity)
    xmin = float(bb.MinPoint.X)
    xmax = float(bb.MaxPoint.X)
    ymin = float(bb.MinPoint.Y)
    ymax = float(bb.MaxPoint.Y)
    return xmin, xmax, ymin, ymax

def _bbox_contains_xy(face, x, y, tol=1e-6):
    xmin, xmax, ymin, ymax = _face_bbox_xy(face)
    return (xmin - tol <= x <= xmax + tol) and (ymin - tol <= y <= ymax + tol)

def _safe_box_xy_info(box):
    """
    Robust XY bbox extraction across SpaceClaim builds.
    Returns (cx, cy, xmin, xmax, ymin, ymax) or raises.
    """
    cx = float(box.Center.X)
    cy = float(box.Center.Y)

    # Try common bounding-box APIs
    if hasattr(box, "MaxPoint") and hasattr(box, "MinPoint"):
        xmin = float(box.MinPoint.X)
        xmax = float(box.MaxPoint.X)
        ymin = float(box.MinPoint.Y)
        ymax = float(box.MaxPoint.Y)
        return cx, cy, xmin, xmax, ymin, ymax

    if hasattr(box, "MaxCorner") and hasattr(box, "MinCorner"):
        xmin = float(box.MinCorner.X)
        xmax = float(box.MaxCorner.X)
        ymin = float(box.MinCorner.Y)
        ymax = float(box.MaxCorner.Y)
        return cx, cy, xmin, xmax, ymin, ymax

    if hasattr(box, "Size"):
        sx = float(box.Size.X) * 0.5
        sy = float(box.Size.Y) * 0.5
        xmin = cx - sx
        xmax = cx + sx
        ymin = cy - sy
        ymax = cy + sy
        return cx, cy, xmin, xmax, ymin, ymax

    if hasattr(box, "Extents"):
        sx = float(box.Extents.X)
        sy = float(box.Extents.Y)
        xmin = cx - sx
        xmax = cx + sx
        ymin = cy - sy
        ymax = cy + sy
        return cx, cy, xmin, xmax, ymin, ymax

    raise Exception("Unsupported Box API: no MaxPoint/MinPoint, MaxCorner/MinCorner, Size, or Extents")

def extrude_from_explicit_curves(name, length_mm, curves):
    root = GetRootPart()
    before = list(root.Bodies)

    if not curves:
        raise Exception("No curves provided for '%s'." % name)

    Fill.Execute(CurveSelection.Create(curves))

    after = list(root.Bodies)
    new_bodies = [b for b in after if b not in before]

    if not new_bodies:
        raise Exception("Explicit Fill created no body for '%s'." % name)

    # pick the intended face from these sheet bodies
    body = new_bodies[-1]
    faces = list(body.Faces)
    if not faces:
        raise Exception("Filled body for '%s' has no faces." % name)

    face = max(faces, key=lambda f: f.Area)   # replace later if needed

    opts = ExtrudeFaceOptions()
    opts.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirZ, MM(length_mm), opts)

    created = list(res.CreatedBodies)
    if not created:
        raise Exception("Extrude produced no result for '%s'." % name)

    new_body = created[0]
    new_body.SetName(name)

    # IMPORTANT: delete temporary sheet-fill bodies
    try:
        Delete.Execute(BodySelection.Create(new_bodies))
    except:
        pass

    # OPTIONAL: clear sketch curves too
    try:
        clear_all_sketch_curves()
    except:
        pass

    # IMPORTANT: return to solid/display mode
    try:
        ViewHelper.SetViewMode(InteractionMode.Solid)
    except:
        pass

    return new_body

def find_body_by_name(name):
    root = GetRootPart()
    for b in list(root.Bodies):
        if getattr(b, "Name", "") == name:
            return b
    for comp in list(root.Components):
        try:
            for b in list(comp.Content.Bodies):
                if getattr(b, "Name", "") == name:
                    return b
        except:
            pass
    return None
# -----------------------------

def create_ns(name, body_or_list):
    """Creates a Named Selection in the Groups tab for ANSYS Mechanical (idempotent)."""

    if isinstance(body_or_list, list):
        items = body_or_list
    else:
        items = [body_or_list]

    # Avoid empty groups
    items = [b for b in items if b is not None]
    if not items:
        print("create_ns skipped (no bodies) for:", name)
        return

    ns_name = "NS_" + name

    # IMPORTANT: delete any existing NS with same name first
    _delete_named_selection_if_exists(ns_name)

    sel = BodySelection.Create(items)
    res = NamedSelection.Create(sel, Selection.Empty())

    # Robust rename
    try:
        res.CreatedNamedSelection.SetName(ns_name)
    except Exception as e:
        print("WARNING: failed to rename Named Selection to", ns_name, ":", e)
        # try to delete the stray "Group1" we just made
        try:
            Delete.Execute(Selection.Create(res.CreatedNamedSelection))
        except:
            pass
        
################# For webcutting ########################################
def make_yz_sheet_face(x_cut=0.0, name="TMP_YZ_SHEET"):
    root = GetRootPart()

    # Big enough to cover everything
    W_ow = W_outer + 2.0*t_overwrap + (2.0*D_drain if drain_opt else 0.0)
    H_ow = H_outer + 2.0*t_overwrap
    halfY = 2.0 * max(W_ow, H_ow)
    halfZ = 2.0 * L_extrude

    # Sketch plane = YZ @ X=x_cut
    frame = Frame.Create(
        Point.Create(MM(x_cut), MM(0), MM(0)),
        Direction.DirY,
        Direction.DirZ
    )
    ViewHelper.SetSketchPlane(Plane.Create(frame))

    # Rectangle in (Y,Z) sketch coords
    SketchLine.Create(P2(-halfY, -halfZ), P2(+halfY, -halfZ))
    SketchLine.Create(P2(+halfY, -halfZ), P2(+halfY, +halfZ))
    SketchLine.Create(P2(+halfY, +halfZ), P2(-halfY, +halfZ))
    SketchLine.Create(P2(-halfY, +halfZ), P2(-halfY, -halfZ))

    before = list(root.Bodies)

    # Prefer Fill (creates planar body from closed sketch)
    try:
        Fill.Execute(Selection.Empty())
    except:
        ViewHelper.SetViewMode(InteractionMode.Solid)

    after = list(root.Bodies)
    new_bodies = [b for b in after if b not in before]
    if not new_bodies:
        raise Exception("make_yz_sheet_face: Fill produced no body (sketch not closed?)")

    sheet = new_bodies[-1]
    sheet.SetName(name)

    # Take its largest face (the planar rectangle)
    faces = list(sheet.Faces)
    if not faces:
        raise Exception("make_yz_sheet_face: sheet has no faces")
    face = max(faces, key=lambda f: f.Area)

    return sheet, face

def split_by_yz_plane_face(body, x_cut=0.0, left_name="LEFT", right_name="RIGHT"):
    if body is None:
        raise Exception("split_by_yz_plane_face: body is None")

    sheet, face = make_yz_sheet_face(x_cut=x_cut)

    # Prefix so we can find pieces even if SC renames stuff
    base = body.Name if getattr(body, "Name", None) else "Body"
    prefix = base + "__SPLIT__"
    body.SetName(prefix + "0")

    # Use the SINGLE face as cutter (no thickness => no gap)
    SplitBody.ByCutter(Selection.Create(body), Selection.Create([face]), True)

    # Delete cutter sheet right away
    try:
        Delete.Execute(Selection.Create(sheet))
    except:
        pass

    # Collect pieces by prefix (root + components)
    pieces = []
    root = GetRootPart()
    for b in list(root.Bodies):
        if getattr(b, "Name", "").startswith(prefix):
            pieces.append(b)
    for comp in list(root.Components):
        try:
            for b in list(comp.Content.Bodies):
                if getattr(b, "Name", "").startswith(prefix):
                    pieces.append(b)
        except:
            pass

    if len(pieces) < 2:
        raise Exception("Split failed: did not find 2+ split pieces.")

    # Sort by X-center: left first, right last
    pieces.sort(key=lambda b: float(b.Shape.GetBoundingBox(Matrix.Identity).Center.X))
    left_piece, right_piece = pieces[0], pieces[-1]

    # Delete any extra fragments if they exist
    extras = [b for b in pieces if (b is not left_piece and b is not right_piece)]
    if extras:
        try:
            Delete.Execute(BodySelection.Create(extras))
        except:
            pass

    left_piece.SetName(left_name)
    right_piece.SetName(right_name)
    return left_piece, right_piece
######################################################################

# 3. BUILD GEOMETRY

# (1) Conductors


c1 = find_body_by_name("conductor[1]")
if c1 is None:
    set_sketch_plane_xy()
    clear_all_sketch_curves()   # see below
    sketch_circle(-dx_cond, 0, r_cond)
    c1 = extrude_and_name("conductor[1]", L_extrude)


c2 = find_body_by_name("conductor[2]")
if c2 is None:
    set_sketch_plane_xy()
    clear_all_sketch_curves()   # see below
    sketch_circle(dx_cond, 0, r_cond)
    c2 = extrude_and_name("conductor[2]", L_extrude)

# (2) Cores (Insulation)

if not (is_a_doublet):

    if (is_a_multilumen):

        w = multilumen_wall_thickness
        total_gap = r_core - r_cond

        # Pitch Circle Radius (centered in the insulation)

        rp = r_cond + (total_gap / 2.0)

        # 1. Calculate Hole Diameter (Dh) based on Angular Spacing
        # The distance between centers along the arc must be approx (Dh + w)
        # Chord length formula: Chord = 2 * rp * sin(pi / n)
        # We want: Chord - Dh = w  =>  Dh = Chord - w

        if n_multilumen_cavity > 1:
            chord = 2.0 * rp * math.sin(math.pi / n_multilumen_cavity)
            dh_angular = chord - w

        else:

            dh_angular = total_gap - 2.0 * w # Fallback for single hole

        # 2. Calculate Hole Diameter (Dh) based on Radial Spacing
        dh_radial = total_gap - 2.0 * w

        # Final Diameter must be the smaller of the two to respect all walls
        d_hole = min(dh_angular, dh_radial)
        r_hole = d_hole / 2.0


        def draw_mirrored_core(center_x, is_right_side):

            set_sketch_plane_xy()
            sketch_circle(center_x, 0, r_core)
            sketch_circle(center_x, 0, r_cond)

            if n_multilumen_cavity > 0 and d_hole > 0:

                for i in range(n_multilumen_cavity):
                    # Calculate angle (offset by 90 deg or pi/n to avoid C2C axis if desired)
                    angle = (2.0 * math.pi * i) / n_multilumen_cavity

                    # Enforce Symmetry
                    if not is_right_side:
                        hx = center_x - rp * math.cos(angle)
                    else:
                        hx = center_x + rp * math.cos(angle)

                    hy = rp * math.sin(angle)
                    sketch_circle(hx, hy, r_hole)

       # Build Cores
       # draw_mirrored_core(-dx_cond, False)
       # extrude_and_name("single_core[1]", L_extrude, True)
       # draw_mirrored_core(dx_cond, True)
       # extrude_and_name("single_core[2]", L_extrude, True)


        # Function to draw the core with lumens
        def merge_bodies_keep_name(name, bodies):
            bodies = [b for b in bodies if b is not None]
            if len(bodies) < 2:
                return bodies[0] if bodies else None

            res = Combine.Merge(BodySelection.Create(bodies))

            # Return shape varies by build
            out = None
            if hasattr(res, "MergedBody") and res.MergedBody is not None:
                out = res.MergedBody
            else:
                try:
                    out = list(res.CreatedBodies)[0]
                except:
                    out = bodies[0]

            try:
                out.SetName(name)
            except:
                pass

            return out

        def sketch_core_with_lumens(center_x):
            set_sketch_plane_xy()
            # 1. Draw the main boundaries of the insulation
            sketch_circle(center_x, 0, r_core) 
            sketch_circle(center_x, 0, r_cond)
            
            if n_multilumen_cavity > 0:
                t_septum = multilumen_wall_thickness
                is_right = (center_x > 0)
                
                # Radial limits for uniform top/bottom walls
                ri = r_cond + t_septum
                ro = r_core - t_septum
                
                angle_step = (2.0 * math.pi) / n_multilumen_cavity
                # The 'mirror' alignment you requested
                base_angle = math.pi if is_right else 0.0

                for i in range(n_multilumen_cavity):
                    # Center axis of the current lumen
                    phi = base_angle + (i * angle_step)
                    
                    if multilumen_shape_opt == 0: # CIRCULAR
                        rp = (ri + ro) / 2.0
                        hx, hy = center_x + rp * math.cos(phi), rp * math.sin(phi)
                        # Ensure circular hole fits in the annulus and the chord space
                        chord = 2.0 * rp * math.sin(angle_step / 2.0)
                        d_hole = min(ro - ri, chord - t_septum)
                        sketch_circle(hx, hy, d_hole / 2.0)
                        
                    else: # PARALLEL-WALLED LUMIN (Trapezoidal-ish)
                        # We need the angle alpha that subtends half the septum thickness at a given radius
                        # To keep walls parallel, we use a fixed distance (t_septum/2) from the center line
                        d_half = t_septum / 2.0
                        
                        # Calculate angular positions at the inner and outer radius to keep sides parallel
                        # theta = asin(distance_from_centerline / radius)
                        alpha_i = math.asin(d_half / ri)
                        alpha_o = math.asin(d_half / ro)
                        
                        # Total angle per segment is angle_step. 
                        # To keep a web of t_septum, the cavity walls are at:
                        # (angle_step/2) - (angular_offset_of_half_septum)
                        
                        # Corners relative to the center_x
                        # Right wall of lumen
                        p1_ang = phi + (angle_step/2.0 - alpha_i)
                        p4_ang = phi + (angle_step/2.0 - alpha_o)
                        # Left wall of lumen
                        p2_ang = phi - (angle_step/2.0 - alpha_i)
                        p3_ang = phi - (angle_step/2.0 - alpha_o)

                        # Convert to XY
                        p1 = P2(center_x + ri * math.cos(p1_ang), ri * math.sin(p1_ang))
                        p2 = P2(center_x + ri * math.cos(p2_ang), ri * math.sin(p2_ang))
                        p3 = P2(center_x + ro * math.cos(p3_ang), ro * math.sin(p3_ang))
                        p4 = P2(center_x + ro * math.cos(p4_ang), ro * math.sin(p4_ang))

                        # Draw the closed loop for the lumen
                        # Arcs must follow the circle (P2(center_x, 0) is the local center)
                        SketchArc.CreateSweepArc(P2(center_x, 0), p1, p2, True)  # Inner Arc
                        SketchLine.Create(p2, p3)                              # Side Wall 1
                        SketchArc.CreateSweepArc(P2(center_x, 0), p3, p4, False) # Outer Arc
                        SketchLine.Create(p4, p1)                              # Side Wall 2

        # Build Core 1
        if core_mode == "separate":
            clear_all_sketch_curves()  
            sketch_core_with_lumens(-dx_cond)
            score1 = extrude_and_name("single_core[1]", L_extrude, True)

            clear_all_sketch_curves() 
            sketch_core_with_lumens(dx_cond)
            score2 = extrude_and_name("single_core[2]", L_extrude, True)


        elif core_mode == "merged":
            # Build BOTH multilumen cores first (left + right), then Boolean-union them
            core1 = find_body_by_name("single_core[1]")
            core2 = find_body_by_name("single_core[2]")
            if (core1 is None) or (core2 is None):
                set_sketch_plane_xy()
                clear_all_sketch_curves() 
                sketch_core_with_lumens(-dx_cond)
                core1 = extrude_and_name("single_core[1]", L_extrude, True)
                set_sketch_plane_xy()
                clear_all_sketch_curves()
                sketch_core_with_lumens(+dx_cond)
                core2 = extrude_and_name("single_core[2]", L_extrude, True)
                #score = extrude_and_name("single_core_merged", L_extrude, True)
                # Union
                
                score = find_body_by_name("single_core_merged")
                if score is None:
                    clear_all_sketch_curves()
                    score = union_bodies("single_core_merged", [core1, core2])


    else:

        if core_mode == "separate":
            set_sketch_plane_xy()
            sketch_circle(-dx_cond, 0, r_core)
            sketch_circle(-dx_cond, 0, r_cond)
            score1 = extrude_and_name("single_core[1]", L_extrude, True)


            set_sketch_plane_xy()
            sketch_circle(dx_cond, 0, r_core)
            sketch_circle(dx_cond, 0, r_cond)
            score2 = extrude_and_name("single_core[2]", L_extrude, True)


        elif core_mode == "merged":
            # Create both overlapping single core bodies first
            score = find_body_by_name("single_core_merged")
            if score is None:
                core1 = find_body_by_name("single_core[1]")
                core2 = find_body_by_name("single_core[2]")

                if core1 is None:
                    set_sketch_plane_xy()
                    clear_all_sketch_curves()
                    sketch_circle(-dx_cond, 0, r_core)
                    sketch_circle(-dx_cond, 0, r_cond)
                    core1 = extrude_and_name("single_core[1]", L_extrude, True)

                if core2 is None:
                    set_sketch_plane_xy()
                    clear_all_sketch_curves()
                    sketch_circle(dx_cond, 0, r_core)
                    sketch_circle(dx_cond, 0, r_cond)
                    core2 = extrude_and_name("single_core[2]", L_extrude, True)

                clear_all_sketch_curves()
                score = union_bodies("single_core_merged", [core1, core2])

                # optional cleanup of originals
                victims = []
                for b in [core1, core2]:
                    if b is not None and b is not score:
                        victims.append(b)
                if victims:
                    try:
                        Delete.Execute(BodySelection.Create(victims))
                    except:
                        pass
            

# (3) Filler / Doublet

set_sketch_plane_xy()

if is_a_doublet:

    sketch_profile(dx_in, R_in, W_in, H_in)
    sketch_circle(-dx_cond, 0, r_cond)
    sketch_circle(dx_cond, 0, r_cond)
    second_extrusion = extrude_and_name("Second_Extrusion", L_extrude, True)
    #create_ns("Second_Extrusion", second_extrusion)

elif filler_mode == "fill":
    
    second_extrusion = find_body_by_name("Second_Extrusion")
    if second_extrusion is None:
        clear_all_sketch_curves() 
        sketch_profile(dx_in, R_in, W_in, H_in)
        sketch_circle(-dx_cond, 0, r_core)
        sketch_circle(dx_cond, 0, r_core)
        second_extrusion = extrude_and_name("Second_Extrusion", L_extrude, True)
    #create_ns("Second_Extrusion", second_extrusion)

elif filler_mode == "shell":

    t_filler_shell = (dx_in + R_in) - (dx_cond + r_core)
    sketch_profile(dx_in, R_in, W_in, H_in)
    W_fill_inner = W_in - 2.0 * t_filler_shell
    H_fill_inner = H_in - 2.0 * t_filler_shell
    R_fill_inner = R_in - t_filler_shell

    if H_fill_inner > 0:
        sketch_profile(dx_in, R_fill_inner, W_fill_inner, H_fill_inner)

    second_extrusion = extrude_and_name("Second_Extrusion", L_extrude, True)
    #create_ns("Second_Extrusion", second_extrusion)

# (4) Shield
print("Shield params: W_outer=%.4f H_outer=%.4f dx_shield=%.4f R_shield=%.4f"
      % (W_outer, H_outer, dx_shield, R_shield))

shield = find_body_by_name("Shield")
shield_outer = find_body_by_name("Shield_outer_tmp")
shield_inner = find_body_by_name("Shield_inner_tmp")


if shield is None:
    if shield_outer is None:
        set_sketch_plane_xy()
        clear_all_sketch_curves()
        outer_curves = sketch_profile(dx_shield, R_shield, W_outer, H_outer)
        shield_outer = extrude_from_explicit_curves("Shield_outer_tmp", L_extrude, outer_curves)
    if shield_inner is None:
        set_sketch_plane_xy()
        clear_all_sketch_curves()
        inner_curves = sketch_profile(dx_in, R_in, W_in, H_in)
        shield_inner = extrude_from_explicit_curves("Shield_inner_tmp", L_extrude, inner_curves)

    options = MakeSolidsOptions()
    options.SubtractFromTarget = True

    Combine.Intersect(
        BodySelection.Create([shield_outer]),
        BodySelection.Create([shield_inner]),
        options
    )

    try:
        Delete.Execute(BodySelection.Create([shield_inner]))
    except:
        pass

    shield = find_body_by_name("Shield_outer_tmp")
    if shield is None:
        shield = shield_outer

    try:
        shield.SetName("Shield")
    except:
        pass


# (5) and (6) Drains & Overwrap
x_drain = dx_shield + R_shield + r_drain
if drain_opt:
    drain1 = find_body_by_name("drain[1]")
    
    if drain1 is None:
        set_sketch_plane_xy()
        clear_all_sketch_curves()
        sketch_circle(-x_drain, 0, r_drain)
        drain1 = extrude_and_name("drain[1]", L_extrude)
        #create_ns("drain[1]", drain1)
    drain2 = find_body_by_name("drain[2]")
    if drain2 is None:
        set_sketch_plane_xy()
        clear_all_sketch_curves()
        sketch_circle(x_drain, 0, r_drain)
        drain2 = extrude_and_name("drain[2]", L_extrude)
        #create_ns("drain[2]", drain2)
        
    # Existing drain-wrap logic relies on Double-D tangency; 
    # This remains as previously defined
    #set_sketch_plane_xy()
    #sketch_wrap_loop_with_drains(0)
    #sketch_wrap_loop_with_drains(t_overwrap)
    #overwrap = extrude_and_name("Overwrap", L_extrude, False)
    overwrap = find_body_by_name("Overwrap")
    if overwrap is None:
        overwrap_outer = find_body_by_name("Overwrap_outer_tmp")
        
        if overwrap_outer is None:
            set_sketch_plane_xy()
            clear_all_sketch_curves()
            outer_curves = sketch_wrap_loop_with_drains(t_overwrap)
            overwrap_outer = extrude_from_explicit_curves("Overwrap_outer_tmp", L_extrude, outer_curves)
        
        overwrap_inner = find_body_by_name("Overwrap_inner_tmp")
        if overwrap_inner is None:
            set_sketch_plane_xy()
            clear_all_sketch_curves()
            inner_curves = sketch_wrap_loop_with_drains(0)
            overwrap_inner = extrude_from_explicit_curves("Overwrap_inner_tmp", L_extrude, inner_curves)

        options = MakeSolidsOptions()
        options.SubtractFromTarget = True

        Combine.Intersect(
            BodySelection.Create([overwrap_outer]),
            BodySelection.Create([overwrap_inner]),
            options
        )

        try:
            live_inner = find_body_by_name("Overwrap_inner_tmp")
            if live_inner is not None:
                Delete.Execute(BodySelection.Create([live_inner]))
        except:
            pass

        overwrap = find_body_by_name("Overwrap_outer_tmp")
        if overwrap is None:
            overwrap = overwrap_outer

        if overwrap is None:
            raise Exception("Could not refetch final Overwrap body after boolean.")

        overwrap.SetName("Overwrap")
else:
    # 1. Create the Outer Solid for Overwrap
    overwrap = find_body_by_name("Overwrap")
    overwrap_outer = find_body_by_name("Overwrap_outer_tmp")
    overwrap_inner = find_body_by_name("Overwrap_inner_tmp")
    if overwrap is None:
        if overwrap_outer is None:
            set_sketch_plane_xy()
            clear_all_sketch_curves()
            outer_curves = sketch_profile(dx_shield, R_shield + t_overwrap, W_outer + 2*t_overwrap, H_outer + 2*t_overwrap)
            overwrap_outer = extrude_from_explicit_curves("Overwrap_outer_tmp", L_extrude, outer_curves)
        if overwrap_inner is None:
            set_sketch_plane_xy()
            clear_all_sketch_curves()
            inner_curves = sketch_profile(dx_shield, R_shield, W_outer, H_outer)
            overwrap_inner = extrude_from_explicit_curves("Overwrap_inner_tmp", L_extrude, inner_curves)

        options = MakeSolidsOptions()
        options.SubtractFromTarget = True

        Combine.Intersect(
            BodySelection.Create([overwrap_outer]),
            BodySelection.Create([overwrap_inner]),
            options
        )

        try:
            Delete.Execute(BodySelection.Create([overwrap_inner]))
        except:
            pass

        overwrap = find_body_by_name("Overwrap_outer_tmp")
        if overwrap is None:
            overwrap = overwrap_outer

        try:
            overwrap.SetName("Overwrap")
        except:
            pass
# -----------------------------

# 4. CLEANUP (Simple Name-based)

# ----------------------------
def cleanup_sheet_bodies(delete_all_non_solids=True, name_filter=None, verbose=True):
    root = GetRootPart()
    sheet_bodies = []

    for b in list(root.Bodies):
        try:
            is_solid = getattr(b, "IsSolid", None)
            # Some versions expose IsSheetBody; use it if available
            is_sheet = getattr(b, "IsSheetBody", None)

            if is_sheet is None:
                # fallback: treat "not solid" as sheet/surface
                is_sheet = (is_solid is False)

            if delete_all_non_solids:
                candidate = (is_sheet is True) or (is_solid is False)
            else:
                candidate = (is_solid is False)

            if candidate:
                if (name_filter is None) or (b.Name == name_filter):
                    sheet_bodies.append(b)

                    if verbose:
                        print("FOUND non-solid body:",
                              "Name='{}'".format(b.Name),
                              "IsSolid={}".format(is_solid),
                              "IsSheetBody={}".format(getattr(b, "IsSheetBody", "NA")))

        except Exception as e:
            if verbose:
                print("Skipping body due to exception:", e)

    if not sheet_bodies:
        if verbose:
            print("No sheet/surface bodies found to delete.")
        return

    # IMPORTANT: Delete expects a proper selection (BodySelection works reliably)
    sel = BodySelection.Create(sheet_bodies)
    Delete.Execute(sel)

    if verbose:
        print("Deleted {} sheet/surface bodies.".format(len(sheet_bodies)))


# Call at very end:
cleanup_sheet_bodies(delete_all_non_solids=True, name_filter=None, verbose=True)

# -----------------------------
# ORGANIZING: Move Cable to Component
# -----------------------------
def _comp_by_name(name):
    root = GetRootPart()
    for c in list(root.Components):
        try:
            if c.GetName() == name:
                return c
        except:
            if getattr(c, "Name", "") == name:
                return c
    return None

def move_all_root_bodies_to_component(comp_name):

    # 1. Get all bodies currently in the Root
    # We convert to a Python list immediately so we don't iterate over a changing collection
    root_bodies = list(GetRootPart().Bodies)
    
    if not root_bodies:
        print("No bodies in Root to move.")
        return
    # 2. Create the new Component inside the Root
    # This creates a new internal component
     # 1. Get the current document
    doc = Window.ActiveWindow.Document

    # 2. Create the Part Definition (Fixing your error)
    # This creates the "blueprint" for the part, but doesn't show it in the tree yet
    new_part_definition = Part.Create(doc, comp_name)

    # 3. Create the Component (Instance)
    # This actually adds it to the assembly under the Root Part
    new_comp = Component.Create(GetRootPart(), new_part_definition)
    
    # 3. Move the bodies
    # ComponentHelper.MoveBodiesToComponent takes (BodySelection, TargetComponent)
    sel = BodySelection.Create(root_bodies)
    ComponentHelper.MoveBodiesToComponent(sel, new_comp)
    
    print("Moved " + str(len(root_bodies)) + " bodies to component: " + comp_name)


# Delete surface bodies by near-zero volume
def delete_bodies_by_exact_name_in_component(comp_name, exact_name="Surface", verbose=True):
    comp = _comp_by_name(comp_name)
    if comp is None:
        print("Component not found:", comp_name)
        return

    victims = []
    for b in list(comp.Content.Bodies):
        nm = getattr(b, "Name", "")
        if nm == exact_name:
            victims.append(b)

    if not victims:
        if verbose:
            print("No bodies named '%s' found in %s." % (exact_name, comp_name))
        return

    Delete.Execute(BodySelection.Create(victims))
    if verbose:
        print("Deleted %d bodies named '%s' in %s." % (len(victims), exact_name, comp_name))


# --- EXECUTE THE MOVE ---
move_all_root_bodies_to_component("cable_bodies")
delete_bodies_by_exact_name_in_component("cable_bodies", "Surface", verbose=True)
    
# ============================================================
# 3-POINT BENDING RIG (GOOD + BAD) — SINGLE SCRIPT
# Choose mode via bending_mode = "good_3point" or "bad_3point" or "good_2pullyes", or "bad_2pulleys"
# ============================================================
# ============================================================
# Helpers (stable / minimal)
# ============================================================

def _body_name(b):
    try: return b.GetName()
    except: return getattr(b, "Name", "")

def get_or_create_component(name):
    root = GetRootPart()
    for c in list(root.Components):
        try:
            if c.GetName() == name: return c
        except:
            if getattr(c, "Name", "") == name: return c

    doc = Window.ActiveWindow.Document
    part_def = Part.Create(doc, name)
    comp = Component.Create(root, part_def)
    return comp

def move_body_to_component(body, target_component):
    nm = _body_name(body)
    sel = BodySelection.Create([body])
    ComponentHelper.MoveBodiesToComponent(sel, target_component)

    # Refetch by name in target
    for b in list(target_component.Content.Bodies):
        if _body_name(b) == nm:
            return b
    # fallback
    return target_component.Content.Bodies[-1]


def move_body_to_component_once(body, current_comp, target_comp, final_name):
    if current_comp is target_comp:
        try:
            body.SetName(final_name)
        except:
            pass
        return body

    tmp_name = final_name + "__TEMP__"
    body.SetName(tmp_name)

    sel = BodySelection.Create([body])
    ComponentHelper.MoveBodiesToComponent(sel, target_comp)

    # Refetch + rename
    for b in list(target_comp.Content.Bodies):
        if _body_name(b) == tmp_name:
            b.SetName(final_name)
            return b

    # Fallback: if it already exists somehow
    for b in list(target_comp.Content.Bodies):
        if _body_name(b) == final_name:
            return b

    raise Exception("Move succeeded but could not refetch moved body: " + final_name)

def find_body_anywhere_by_name(name):
    root = GetRootPart()
    for b in list(root.Bodies):
        if _body_name(b) == name:
            return b, None
    for comp in list(root.Components):
        try:
            for b in list(comp.Content.Bodies):
                if _body_name(b) == name:
                    return b, comp
        except:
            pass
    return None, None

def find_bodies_by_prefix_anywhere(prefix):
    root = GetRootPart()
    out = []
    for b in list(root.Bodies):
        if _body_name(b).startswith(prefix):
            out.append(b)
    for comp in list(root.Components):
        try:
            for b in list(comp.Content.Bodies):
                if _body_name(b).startswith(prefix):
                    out.append(b)
        except:
            pass
    return out

# -----------------------------
# Sketch plane helpers
# -----------------------------
def set_sketch_plane_yz_at_x(x0):
    frame = Frame.Create(
        Point.Create(MM(x0), MM(0), MM(0)),
        Direction.DirY,
        Direction.DirZ
    )
    ViewHelper.SetSketchPlane(Plane.Create(frame))

def set_sketch_plane_xy_at_z(z0):
    frame = Frame.Create(
        Point.Create(MM(0), MM(0), MM(z0)),
        Direction.DirX,
        Direction.DirY
    )
    ViewHelper.SetSketchPlane(Plane.Create(frame))

def set_sketch_plane_zx_at_y(y0):
    # ZX plane at Y=y0. In-plane axes: Z (u), X (v)
    frame = Frame.Create(
        Point.Create(MM(0), MM(y0), MM(0)),
        Direction.DirZ,
        Direction.DirX
    )
    ViewHelper.SetSketchPlane(Plane.Create(frame))

# -----------------------------
# Sketch primitives
# -----------------------------
def sketch_rectangle_xy(xmin, xmax, ymin, ymax):
    SketchLine.Create(P2(xmin, ymin), P2(xmax, ymin))
    SketchLine.Create(P2(xmax, ymin), P2(xmax, ymax))
    SketchLine.Create(P2(xmax, ymax), P2(xmin, ymax))
    SketchLine.Create(P2(xmin, ymax), P2(xmin, ymin))

def sketch_rectangle_uv(umin, umax, vmin, vmax):
    # Rectangle in current sketch plane coords (u,v) mapped to Point2D(x,y)
    SketchLine.Create(P2(umin, vmin), P2(umax, vmin))
    SketchLine.Create(P2(umax, vmin), P2(umax, vmax))
    SketchLine.Create(P2(umax, vmax), P2(umin, vmax))
    SketchLine.Create(P2(umin, vmax), P2(umin, vmin))

def sketch_doubleD_shifted(x0, y0, dx, R):
    SketchArc.CreateSweepArc(P2(x0 + dx, y0), P2(x0 + dx, y0 + R), P2(x0 + dx, y0 - R), True)
    SketchArc.CreateSweepArc(P2(x0 - dx, y0), P2(x0 - dx, y0 + R), P2(x0 - dx, y0 - R), False)
    SketchLine.Create(P2(x0 - dx, y0 + R), P2(x0 + dx, y0 + R))
    SketchLine.Create(P2(x0 - dx, y0 - R), P2(x0 + dx, y0 - R))

# -----------------------------
# Extrude helpers
# -----------------------------
def extrude_last_profile_along_x(name, length_mm, pick_largest=True):
    solidify_sketch()
    temp_body = GetRootPart().Bodies[-1]
    face = max(temp_body.Faces, key=lambda f: f.Area) if pick_largest else temp_body.Faces[0]
    opts = ExtrudeFaceOptions()
    opts.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirX, MM(length_mm), opts)
    new_body = list(res.CreatedBodies)[0]
    new_body.SetName(name)
    return new_body

def extrude_last_profile_along_y(name, length_mm, pick_largest=True):
    solidify_sketch()
    temp_body = GetRootPart().Bodies[-1]
    face = max(temp_body.Faces, key=lambda f: f.Area) if pick_largest else temp_body.Faces[0]
    opts = ExtrudeFaceOptions()
    opts.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirY, MM(length_mm), opts)
    new_body = list(res.CreatedBodies)[0]
    new_body.SetName(name)
    return new_body

def extrude_largest_face_along_z(name, length_mm):
    solidify_sketch()
    temp_body = GetRootPart().Bodies[-1]
    faces = list(temp_body.Faces)
    if not faces:
        raise Exception("No faces found to extrude for %s (sketch not closed?)" % name)
    face = max(faces, key=lambda f: f.Area)
    opts = ExtrudeFaceOptions()
    opts.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirZ, MM(length_mm), opts)
    created = list(res.CreatedBodies)
    if not created:
        raise Exception("Extrude created no bodies for %s" % name)
    body = created[0]
    body.SetName(name)
    return body



# --- NOW create Named Selections (last) ---
create_ns("conductor[1]", find_body_anywhere_by_name("conductor[1]")[0])
create_ns("conductor[2]", find_body_anywhere_by_name("conductor[2]")[0])
if core_mode == "separate":
    create_ns("single_core[1]", find_body_anywhere_by_name("single_core[1]")[0])
    create_ns("single_core[2]", find_body_anywhere_by_name("single_core[2]")[0])
else:
    create_ns("single_core_merged", find_body_anywhere_by_name("single_core_merged")[0])
create_ns("Second_Extrusion", find_body_anywhere_by_name("Second_Extrusion")[0])
create_ns("Shield", find_body_anywhere_by_name("Shield")[0])
if drain_opt:
    create_ns("drain[1]", find_body_anywhere_by_name("drain[1]")[0])
    create_ns("drain[2]", find_body_anywhere_by_name("drain[2]")[0])
create_ns("Overwrap", find_body_anywhere_by_name("Overwrap")[0])


# --- SHARE TOPOLOGY / IMPRINT CHAIN (do this BEFORE Named Selections) ---
def _safe_Name(obj):
    # Always return a plain python string, never throw.
    if obj is None:
        return "<None>"

    # 1) Try Name (property) safely
    try:
        n = getattr(obj, "Name", None)
        if n is not None:
            try:
                # sometimes it's already a string, sometimes not
                return str(n)
            except:
                return "<Name-unprintable>"
    except:
        pass

    # 2) Try GetName() safely (some SC objects prefer this)
    try:
        if hasattr(obj, "GetName"):
            n = obj.GetName()
            if n is not None:
                try:
                    return str(n)
                except:
                    return "<GetName-unprintable>"
    except:
        pass

    # 3) Last resort: type info only (avoid str(obj))
    try:
        return "<%s>" % obj.GetType().FullName
    except:
        pass

    return "<UnknownObject>"

def _get_parent_component_or_part(body):
    """
    In SpaceClaim, body.Parent is often a Component (or Part/RootPart-ish depending on context).
    We treat it as the object that may carry SharedTopology.
    """
    try:
        return body.Parent
    except:
        return None

def _set_shared_topology_on_parent(body, verbose=True):
    parent = _get_parent_component_or_part(body)
    if parent is None:
        if verbose:
            print("SharedTopology: could not access parent for", _safe_name(body))
        return False

    # Some builds expose SharedTopology on the parent component; some don’t.
    if hasattr(parent, "SharedTopology"):
        try:
            parent.SharedTopology = SharedTopology.Share  # Share, Merge, None
            if verbose:
                print("SharedTopology set to Share on parent of", _safe_name(body),
                      "parent=", _safe_name(parent))
            return True
        except Exception as e:
            if verbose:
                print("SharedTopology: failed to set on parent of", _safe_name(body), ":", e)
            return False

    if verbose:
        print("SharedTopology: parent has no SharedTopology attribute for", _safe_name(body),
              "parent=", _safe_name(parent))
    return False

def share_topology_pair_WIP(body_a, body_b,
                       mode="share_only",
                       imprint_tol_mm=None,
                       verbose=True):
    """
    mode:
      - "share_only": sets SharedTopology on parents (HFSS-cleanest)
      - "share_and_imprint": sets SharedTopology + tries Imprint/Repair.Share (strongest, but can fragment faces)

    imprint_tol_mm:
      - None: leave default imprint tolerance
      - float: sets ImprintOptions().Tolerance = MM(imprint_tol_mm) when available
    """

    if body_a is None or body_b is None:
        return False

    if body_a is body_b:
        if verbose:
            print("share_topology_pair: same body provided; skipping:", _safe_name(body_a))
        return True

    # 1) Always do SharedTopology flag first (low-risk, Mechanical-friendly)
    ok_a = _set_shared_topology_on_parent(body_a, verbose=verbose)
    ok_b = _set_shared_topology_on_parent(body_b, verbose=verbose)

    # If user only wants the “WB share topology” behavior, stop here.
    if mode.lower().strip() in ["share_only", "share"]:
        if verbose:
            print("share_topology_pair: mode=share_only; done.")
        return (ok_a or ok_b)

    # 2) Optional: physical imprint / repair-share (higher risk for HFSS)
    sel = None
    try:
        sel = BodySelection.Create([body_a, body_b])
    except Exception as e:
        if verbose:
            print("share_topology_pair: failed to create BodySelection:", e)
        return False

    # Try Imprint.Execute first (if available)
    try:
        options = ImprintOptions()
        # Not all builds expose Tolerance on ImprintOptions, so guard it.
        if imprint_tol_mm is not None and hasattr(options, "Tolerance"):
            options.Tolerance = MM(float(imprint_tol_mm))

        Imprint.Execute(sel, options)
        if verbose:
            print("Imprint.Execute OK between:", _safe_name(body_a), "and", _safe_name(body_b),
                  ("tol=%.6g mm" % imprint_tol_mm) if imprint_tol_mm is not None else "")
        return True

    except Exception as e_imprint:
        if verbose:
            print("Imprint.Execute failed:", e_imprint)

    # Fallback: Repair Share tool (sometimes exists as Share.Share)
    try:
        # Depending on the build, Share.Share might accept Selection.Create([...]) or BodySelection.
        try:
            Share.Share(sel)
        except:
            Share.Share(Selection.Create([body_a, body_b]))
        if verbose:
            print("Repair.Share OK between:", _safe_name(body_a), "and", _safe_name(body_b))
        return True
    except Exception as e_share:
        if verbose:
            print("Repair.Share failed:", e_share)

    if verbose:
        print("share_topology_pair: all methods failed for:",
              _safe_name(body_a), "<->", _safe_name(body_b))
    return False

def _as_body(x):
    if x is None: return None
    if isinstance(x, list) or isinstance(x, tuple):
        return x[0] if x else None
    return x    

def share_topology_pair(body_a, body_b, mode="share_only", verbose=True):
    """
    mode:
      - "share_only": sets SharedTopology on parents (HFSS-cleanest)
      - "share_and_imprint": sets SharedTopology + tries Imprint/Repair.Share (strongest, but can fragment faces)
    """
    if body_a is None or body_b is None:
        return

    # 1. Set the 'Share Topology' property on the Component
    # This is what Mechanical 2025 R2 looks at to merge nodes.
    try:
        # Get the component containing body_a
        comp = body_a.Parent
        if hasattr(comp, 'SharedTopology'):
            comp.SharedTopology = SharedTopology.Share # Options: Share, Merge, None
            if verbose: print("Component Property set to 'Share' for:", body_a.Name)
    except Exception as e:
        if verbose: print("Failed to set SharedTopology property:", e)

    if mode != "share_and_imprint":
        return  
    
    # 2. Physical Imprint (The 'Repair' way)
    # This physically creates the edges where bodies touch.
    try:
        selection = BodySelection.Create([body_a, body_b])
        # In 252, Imprint is often found in the 'Repair' namespace or via the 'PowerSelect'
        # This is the most stable 252 command for imprinting
        options = ImprintOptions()
        # You can set options.Tolerance = MM(0.01) if needed
        Imprint.Execute(selection, options)
        
        if verbose: print("Physical Imprint successful between:", body_a.Name, "and", body_b.Name)
    except Exception as e_imprint:
        if verbose: print("Physical Imprint failed:", e_imprint)
        
        # Final Fallback: The 'Share' tool in the Repair tab
        try:
            Share.Share(selection) # This is the 2025 R2 'Share' button logic
            if verbose: print("Repair.Share successful.")
        except:
            if verbose: print("All sharing methods failed.")

def share_topology_group(bodies, do_imprint=False, imprint_tol_mm=None, verbose=True):
    """
    bodies: list of Body objects (e.g. [c1, c2, core])
    do_imprint: if True, run Imprint/Share on the whole set at once
    imprint_tol_mm: optional imprint tolerance in mm (e.g. 0.005)
    """
    bodies = [b for b in bodies if b is not None]
    if len(bodies) < 2:
        return

    # 1) Set SharedTopology (best effort) on the *container* of these bodies.
    # NOTE: your existing code uses body_a.Parent; that's often not the right level.
    # We’ll keep it minimal: try each body's parent.
    try:
        for b in bodies:
            comp = getattr(b, "Parent", None)
            if comp is not None and hasattr(comp, "SharedTopology"):
                comp.SharedTopology = SharedTopology.Share
        if verbose: print("SharedTopology set to Share on parent(s) (best effort).")
    except Exception as e:
        if verbose: print("Failed setting SharedTopology:", e)

    if not do_imprint:
        return

    # 2) Imprint once on ALL bodies
    selection = BodySelection.Create(bodies)
    try:
        opts = ImprintOptions()
        if imprint_tol_mm is not None and hasattr(opts, "Tolerance"):
            opts.Tolerance = MM(imprint_tol_mm)
        Imprint.Execute(selection, opts)
        if verbose: print("Imprint successful on group of", len(bodies), "bodies.")
        return
    except Exception as e_imprint:
        if verbose: print("Imprint failed:", e_imprint)

    # 3) Fallback: Repair > Share once on ALL bodies
    try:
        Share.Share(selection)
        if verbose: print("Repair.Share successful on group.")
    except Exception as e_share:
        if verbose: print("Repair.Share failed:", e_share)

# AFTER all relevant bodies exist and are final (including merged core case)
# Use share_only for HFSS cleanliness

if core_mode != "merged":
    share_topology_pair(c1, score1, mode="share_only", verbose=True)
    share_topology_pair(c2, score2, mode="share_only", verbose=True)
    share_topology_pair(score1, second_extrusion, mode="share_only", verbose=True)
    share_topology_pair(score2, second_extrusion, mode="share_only", verbose=True)
else:
    share_topology_pair(core1, core2, mode="share_only", verbose=True)
    share_topology_pair(c1, core1, mode="share_only", verbose=True)
    share_topology_pair(c2, core2, mode="share_only", verbose=True)
    share_topology_pair(core1, second_extrusion, mode="share_only", verbose=True)
    share_topology_pair(core2, second_extrusion, mode="share_only", verbose=True)


# ---------------- vector helpers ----------------

def _dot(u, v): return u.X*v.X + u.Y*v.Y + u.Z*v.Z
def _norm(u): return math.sqrt(_dot(u,u))
def _scale(u, s): return Direction.Create(u.X*s, u.Y*s, u.Z*s)
def _sub(a, b): return Direction.Create(a.X-b.X, a.Y-b.Y, a.Z-b.Z)

def _unit(u):
    m = _norm(u)
    if m == 0: return Direction.Create(0,0,0)
    return Direction.Create(u.X/m, u.Y/m, u.Z/m)

def _reject(u, axis_unit):
    # remove component along axis
    return _sub(u, _scale(axis_unit, _dot(u, axis_unit)))

def _zero_z(u):  # keep XY only
    return Direction.Create(u.X, u.Y, 0.0)

# ---------------- API-dependent shims ----------------


def _face_point(face):
    """
    Return a representative point on the face using parameter evaluation 
    or bounding box center for SpaceClaim 252.
    """
    try:
        # 1. Attempt to get a point at the middle of the UV domain
        # This is the most 'representative' point on the actual surface
        domain = face.Geometry.Domain
        mid_u = (domain.RangeU.Start + domain.RangeU.End) / 2.0
        mid_v = (domain.RangeV.Start + domain.RangeV.End) / 2.0
        
        # Evaluate returns a FaceEvaluation object containing the Point
        return face.Geometry.Evaluate(Point2D.Create(mid_u, mid_v)).Point
    except:
        try:
            # 2. Fallback: Bounding box center (requires Matrix.Identity in V252)
            # Access the underlying Modeler Face for the BoundingBox
            return face.Shape.GetBoundingBox(Matrix.Identity).Center
        except:
            # 3. Last resort: The position of the first vertex
            return list(face.Vertices)[0].Position

def _face_normal(face, point):
    """
    Evaluates the normal of the face at the point closest to 'point'.
    """
    # 1. Access the Modeler Face through .Shape
    modeler_face = face.Shape 
    
    # 2. Use Geometry.ProjectPoint on the Modeler object
    # This finds the closest point on the surface to your input 'point'
    result = modeler_face.Geometry.ProjectPoint(point)
    
    # 3. Return the Normal vector at that projected location
    return result.Normal

def _get_body_by_exact_name(name):
    """Helper to refetch a body from the root if a pointer dies."""
    for b in GetRootPart().Bodies:
        if b.Name == name: return b
    return None

def is_valid(obj):
    """
    Checks if a SpaceClaim object is still 'alive' in the database.
    """
    if obj is None:
        return False
    # 1. Check the .IsDisposed property (Standard for .NET-based APIs)
    try:
        if obj.IsDisposed:
            return False
    except:
        pass
    
    # 2. Verify it still has a parent (if it's a body, it should be in a Part)
    try:
        if obj.Parent is None:
            return False
    except:
        return False
        
    return True

def get_live_body(body_ref):
    """
    Returns the live version of a body. If the reference is dead,
    it tries to find it by name in the Root Part.
    """
    if is_valid(body_ref):
        return body_ref
    
    # If the pointer is dead, refetch by name
    name_to_find = getattr(body_ref, "Name", None)
    if name_to_find:
        for b in GetRootPart().Bodies:
            if b.Name == name_to_find:
                return b
                
    return None

def body_xy_center_from_bbox(body):
    # Ensure we have a valid reference
    live_body = get_live_body(body)
    if live_body is None:
        raise Exception("Cannot find center: Body reference is dead/deleted.")

    # Access the shape safely
    shape = live_body.Shape
    bbox = shape.GetBoundingBox(Matrix.Identity)
    center_pt = bbox.Center
    return Direction.Create(center_pt.X, center_pt.Y, 0.0)


def _create_named_selection_from_faces(faces, ns_name):
    # Use the logic that worked in your previous scripts
    _delete_named_selection_if_exists(ns_name)
    sel = FaceSelection.Create(faces)
    res = NamedSelection.Create(sel, Selection.Empty())
    if res.CreatedNamedSelection:
        res.CreatedNamedSelection.SetName(ns_name)

# ---------------- main logic ----------------

def create_layer_side_ns_by_normal(layer_body, ns_outer, ns_inner=None,
                                    axis="Z", endcap_cos=0.9, area_min=0.0):
    
    axis_unit = {"X": Direction.DirX, "Y": Direction.DirY, "Z": Direction.DirZ}[axis]
    axis_unit = _unit(axis_unit)

    # Get center for the radial test
    c_xy = body_xy_center_from_bbox(layer_body)
    
    outer_faces = []
    inner_faces = []

    for f in list(layer_body.Faces):
        if area_min > 0.0 and f.Area < area_min:
            continue

        p = _face_point(f)
        n = _face_normal(f, p)
        
        # 1. Exclude end caps (faces parallel to the Z-axis)
        if abs(_dot(n, axis_unit)) > endcap_cos:
            continue

        # 2. Radial test: compare normal vector (n) with position vector (r)
        # r = vector from cable center to face point
        pdir = Direction.Create(p.X, p.Y, p.Z)
        r = _zero_z(_sub(pdir, c_xy))
        
        # If normal points in same direction as radial vector, it's an OUTER face
        if _dot(n, r) >= 0:
            outer_faces.append(f)
        else:
            inner_faces.append(f)

    if outer_faces:
        _create_named_selection_from_faces(outer_faces, ns_outer)
    
    if ns_inner and inner_faces:
        _create_named_selection_from_faces(inner_faces, ns_inner)

    return outer_faces, inner_faces



# -----------------------------
# Parameters (derived defaults)
# -----------------------------
Initial_Gap  = get_param("Initial_Gap",  0.05 * (H_outer + 2.0*t_overwrap))
Nose_Diam    = get_param("Nose_Diam",    1.5  * (H_outer + 2.0*t_overwrap))
Nose_Length  = get_param("Nose_Length",  4.0  * (H_outer + 2.0*t_overwrap))

# ============================================================
# 7A) Loading Nose (GOOD/BAD)
# ============================================================
def create_loading_nose(mode):

    if mode == "good_3point":
        comp = get_or_create_component("RigidParts_3Point_Bending")

        existing, existing_comp = find_body_anywhere_by_name("Rig_Loading_Nose")
        if existing is not None:
            return move_body_to_component_once(existing, existing_comp, comp, "Rig_Loading_Nose")
        
    r = Nose_Diam / 2.0
    z_center = 0.5 * L_extrude

    if mode == "good_3point":
        # Good: nose above (+Y), sketch on YZ@x, extrude along X
        cable_bot_y = - 0.5 * (H_outer + 2.0*t_overwrap)
        y_center = cable_bot_y - Initial_Gap - r

        set_sketch_plane_yz_at_x(-0.5 * Nose_Length)
        sketch_circle(y_center, z_center, r)

        temp_nose = extrude_last_profile_along_x("Rig_Loading_Nose", Nose_Length)
        comp = get_or_create_component("RigidParts_3Point_Bending")
        final_nose = move_body_to_component(temp_nose, comp)
        create_ns("Rig_Loading_Nose", final_nose)
        return final_nose

    elif mode == "bad_3point":
        # Bad: nose on side (+X), sketch on ZX@y, extrude along Y
        cable_bot_x = - 0.5 * (W_outer + 2.0*t_overwrap + (2.0*D_drain if drain_opt else 0.0))
        x_center = cable_bot_x - Initial_Gap - r

        set_sketch_plane_zx_at_y(-0.5 * Nose_Length)
        sketch_circle(z_center, x_center, r)

        temp_nose = extrude_last_profile_along_y("Rig_Loading_Nose", Nose_Length)
        comp = get_or_create_component("RigidParts_3Point_Bending")
        final_nose = move_body_to_component(temp_nose, comp)
        create_ns("Rig_Loading_Nose", final_nose)
        return final_nose

    else:
        raise Exception("Unknown mode: %s (use 'good' or 'bad')" % mode)

# ============================================================
# 7F/7B) Supports (shared geometry; different names/components)
# ============================================================
def cable_envelope_with_optional_drains():
    """
    Returns (W_eff, H_eff) for the OUTER envelope the supports should cover.
    - No drains: classic W_outer/H_outer + 2*t_overwrap
    - With drains: width expands to include drain OD + overwrap
    """
    H_eff = H_outer + 2.0*t_overwrap

    if drain_opt:
        # drain center is already defined in your script before overwrap
        # x_drain = dx_shield + R_shield + r_drain
        half_width = x_drain + (r_drain + t_overwrap)
        W_eff = 2.0 * half_width
    else:
        W_eff = W_outer + 2.0*t_overwrap

    return W_eff, H_eff


def create_support_with_doubleD_pocket(name, z_center):
    W_ow, H_ow = cable_envelope_with_optional_drains()
    R_ow  = 0.5 * H_ow
    dx_ow = 0.5 * (W_ow - H_ow)

    clearance = 0.03 * min(W_ow, H_ow)
    R_pocket  = R_ow + clearance
    dx_pocket = dx_ow

    t_support = 0.1 * L_extrude
    #t_support = max(1.1*H_ow, min(t_support, 2.0*H_ow))

    pad_x    = 0.25 * W_ow
    pad_y_up = pad_x #0.15 * H_ow
    pad_y_dn = pad_x #0.55 * H_ow

    x0 = 0.0
    y0 = 0.0  # your later version

    z0 = z_center - 0.5*t_support
    set_sketch_plane_xy_at_z(z0)

    xmin = -(0.5*W_ow + pad_x)
    xmax = +(0.5*W_ow + pad_x)
    ymin = -(0.5*H_ow + pad_y_dn)
    ymax = +(0.5*H_ow + pad_y_up)
    sketch_rectangle_xy(xmin, xmax, ymin, ymax)

    sketch_doubleD_shifted(x0, y0, dx_pocket, R_pocket)

    return extrude_largest_face_along_z(name, t_support)

def push_support_to_negative_y(body):
    W_ow, H_ow = cable_envelope_with_optional_drains()
    cable_bot_y = -0.5 * H_ow
    # put support top at cable bottom - gap
    target_top_y = cable_bot_y - Initial_Gap

    # current support top y (quick estimate via face sample)
    y_max = None
    for f in list(body.Faces):
        try:
            y = f.Eval(0.5, 0.5).Point.Y
            y_max = y if (y_max is None or y > y_max) else y_max
        except:
            pass
    if y_max is None:
        return body

    dy = target_top_y - y_max
    Move.Translate(Selection.Create(body), Direction.DirY, MM(dy))
    return body

def create_two_supports(mode):
    z1 = 0.10 * L_extrude
    z2 = 0.90 * L_extrude

    if mode in ["good_3point", "bad_3point"]:
        comp_name = "RigidParts_3Point_Bending"
        n1, n2 = "Rig_Support_1", "Rig_Support_2"
    else:
        raise Exception("Unknown mode: %s" % mode)

    rigid_comp = get_or_create_component(comp_name)

    b1, c1 = find_body_anywhere_by_name(n1)
    b2, c2 = find_body_anywhere_by_name(n2)

    if b1 is not None and b2 is not None:
        b1 = move_body_to_component_once(b1, c1, rigid_comp, n1)
        b2 = move_body_to_component_once(b2, c2, rigid_comp, n2)
        return b1, b2

    if b1 is None:
        s1_tmp = create_support_with_doubleD_pocket(n1, z1)
        if mode == "good_3point":
            s1 = split_by_ZX_faces_keep_bottom(s1_tmp, y_cut=0.0, final_name=n1)
        else:
            s1 = split_by_YZ_faces_keep_bottom(s1_tmp, x_cut=0.0, final_name=n1)
        b1 = move_body_to_component_once(s1, None, rigid_comp, n1)

    if b2 is None:
        s2_tmp = create_support_with_doubleD_pocket(n2, z2)
        if mode == "good_3point":
            s2 = split_by_ZX_faces_keep_bottom(s2_tmp, y_cut=0.0, final_name=n2)
        else:
            s2 = split_by_YZ_faces_keep_bottom(s2_tmp, x_cut=0.0, final_name=n2)
        b2 = move_body_to_component_once(s2_tmp, None, rigid_comp, n2)

    return b1, b2
# ============================================================
# 7G) Split support openers (GOOD: ZX@Y=0 keep min Y) (BAD: YZ@X=0 keep min X)
# ============================================================
def _faces_with_normal_near_dir(faces, dir_vec, tol=0.95):
    out = []
    for f in faces:
        try:
            n = f.Plane.Normal
            d = n.X*dir_vec.X + n.Y*dir_vec.Y + n.Z*dir_vec.Z
            if d > tol:
                out.append(f)
        except:
            pass
    return out

def body_representative_y(body):
    # Prefer bounding box when available (stable)
    try:
        bb = body.GetBoundingBox()
        return float(bb.Min.Y)   # smallest Y is what you want for "bottom"
    except:
        pass

    # Fallback to face eval (your old behavior)
    for f in list(body.Faces):
        try:
            return f.Eval(0.5, 0.5).Point.Y
        except:
            pass
    return 0.0

def body_representative_x(body):
    for f in list(body.Faces):
        try:
            return f.Eval(0.5, 0.5).Point.X
        except:
            pass
    return 0.0

def make_thin_cutter_slab_ZX(y_cut=0.0, name="TMP_ZX_CUTTER"):
    W_ow = W_outer + 2.0*t_overwrap
    H_ow = H_outer + 2.0*t_overwrap
    halfZ = 1.20 * L_extrude
    halfX = 1.50 * W_ow
    tY = max(0.02 * H_ow, 0.05)

    set_sketch_plane_zx_at_y(y_cut - 0.5*tY)
    # In ZX sketch: u=Z, v=X
    sketch_rectangle_uv(-halfZ, +halfZ, -halfX, +halfX)

    solidify_sketch()
    temp_body = GetRootPart().Bodies[-1]
    face = max(list(temp_body.Faces), key=lambda f: f.Area)

    opts = ExtrudeFaceOptions()
    opts.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirY, MM(tY), opts)
    slab = list(res.CreatedBodies)[0]
    slab.SetName(name)
    return slab

def make_thin_cutter_slab_YZ(x_cut=0.0, name="TMP_YZ_CUTTER"):
    W_ow = W_outer + 2.0*t_overwrap
    H_ow = H_outer + 2.0*t_overwrap
    halfZ = 1.20 * L_extrude
    halfY = 1.50 * H_ow
    tX = max(0.02 * W_ow, 0.05)

    set_sketch_plane_yz_at_x(x_cut - 0.5*tX)
    # In YZ sketch: u=Z, v=Y (your working convention)
    sketch_rectangle_uv(-halfZ, +halfZ, -halfY, +halfY)

    solidify_sketch()
    temp_body = GetRootPart().Bodies[-1]
    face = max(list(temp_body.Faces), key=lambda f: f.Area)

    opts = ExtrudeFaceOptions()
    opts.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirX, MM(tX), opts)
    slab = list(res.CreatedBodies)[0]
    slab.SetName(name)
    return slab

def split_by_ZX_faces_keep_bottom(body, y_cut=0.0, final_name=None):
    slab = make_thin_cutter_slab_ZX(y_cut=y_cut, name="TMP_ZX_CUTTER")
    slab_faces = list(slab.Faces)

    plusY = _faces_with_normal_near_dir(slab_faces, Direction.DirY)

    try:
        dir_minusY = Direction.DirY.Negate()
    except:
        dir_minusY = Direction.Create(0, -1, 0)

    minusY = _faces_with_normal_near_dir(slab_faces, dir_minusY)

    tool_faces = []
    if plusY:  tool_faces.append(max(plusY,  key=lambda f: f.Area))
    if minusY: tool_faces.append(max(minusY, key=lambda f: f.Area))
    if len(tool_faces) < 2:
        tool_faces = sorted(slab_faces, key=lambda f: f.Area, reverse=True)[:2]

    tag = (final_name if final_name else _body_name(body)) or "RigSupport"
    prefix = tag + "__SPLIT__"
    body.SetName(prefix + "0")

    SplitBody.ByCutter(Selection.Create(body), Selection.Create(tool_faces), True)

    pieces = find_bodies_by_prefix_anywhere(prefix)
    if len(pieces) < 2:
        try: Delete.Execute(Selection.Create(slab))
        except: pass
        if final_name:
            try: body.SetName(final_name)
            except: pass
        print("WARNING: Could not find split pieces by prefix; leaving body as-is.")
        return body

    keeper = None
    keeper_y = None
    for b in pieces:
        yb = body_representative_y(b)
        if keeper is None or yb < keeper_y:
            keeper, keeper_y = b, yb

    Delete.Execute(Selection.Create([b for b in pieces if b is not keeper]))
    Delete.Execute(Selection.Create(slab))

    keeper.SetName(final_name if final_name else tag)
    return keeper

def split_by_YZ_faces_keep_bottom(body, x_cut=0.0, final_name=None):
    slab = make_thin_cutter_slab_YZ(x_cut=x_cut, name="TMP_YZ_CUTTER")
    slab_faces = list(slab.Faces)

    plusX  = _faces_with_normal_near_dir(slab_faces, Direction.DirX)
    try:
        dir_minusX = Direction.DirX.Negate()
    except:
        dir_minusX = Direction.Create(-1, 0, 0)

    minusX = _faces_with_normal_near_dir(slab_faces, dir_minusX)

    tool_faces = []
    if plusX:  tool_faces.append(max(plusX,  key=lambda f: f.Area))
    if minusX: tool_faces.append(max(minusX, key=lambda f: f.Area))
    if len(tool_faces) < 2:
        tool_faces = sorted(slab_faces, key=lambda f: f.Area, reverse=True)[:2]

    tag = (final_name if final_name else _body_name(body)) or "RigSupport"
    prefix = tag + "__SPLIT__"
    body.SetName(prefix + "0")

    SplitBody.ByCutter(Selection.Create(body), Selection.Create(tool_faces), True)

    pieces = find_bodies_by_prefix_anywhere(prefix)
    if len(pieces) < 2:
        try: Delete.Execute(Selection.Create(slab))
        except: pass
        if final_name:
            try: body.SetName(final_name)
            except: pass
        print("WARNING: Could not find split pieces by prefix; leaving body as-is.")
        return body

    keeper = None
    keeper_x = None
    for b in pieces:
        xb = body_representative_x(b)
        if keeper is None or xb < keeper_x:
            keeper, keeper_x = b, xb

    Delete.Execute(Selection.Create([b for b in pieces if b is not keeper]))
    Delete.Execute(Selection.Create(slab))

    keeper.SetName(final_name if final_name else tag)
    return keeper

def angle_360():
    # Try a few common patterns used across SC versions
    try:
        return Angle.Create(2.0 * math.pi)      # radians
    except:
        pass
    try:
        return Angle.Create(MM(360.0))          # unlikely but harmless
    except:
        pass
    try:
        return Angle.Degrees(360.0)             # some builds expose this
    except:
        pass
    try:
        return 360.0                             # some builds just accept float degrees
    except:
        pass
    return 360.0
def pick_planar_face(body):
    # Revolve wants a real planar face; pick the largest is usually correct
    faces = list(body.Faces)
    if not faces:
        raise Exception("No faces found on temp body for revolve (sketch not closed?)")
    return max(faces, key=lambda f: f.Area)
def create_groove_tool_revolved_about_x(name, pulley_R, W_ow, H_ow, clear):
    # Inflate the cable envelope for clearance
    Wg = W_ow + 2.0*clear
    Hg = H_ow + 2.0*clear
    Rg  = 0.5 * Hg
    dxg = 0.5 * (Wg - Hg)

    # Place profile center at (Y=pulley_R, Z=0) in the YZ plane
    yc = pulley_R
    zc = 0.0

    # Sketch on YZ at x = 0
    set_sketch_plane_yz_at_x(0.0)

    # IMPORTANT: your sketch_doubleD_shifted_YZ must actually draw a CLOSED loop
    sketch_doubleD_shifted_YZ(yc, zc, dxg, Rg)

    solidify_sketch()

    # The last created "sketch fill" body is usually at root
    temp = GetRootPart().Bodies[-1]
    face = pick_planar_face(temp)

    # Axis of revolution: X axis through origin
    axis_line = Line.Create(Point.Create(MM(0), MM(0), MM(0)), Direction.DirX)

    opts = RevolveFaceOptions()
    opts.RevolveType = RevolveType.ForceIndependent

    ang = angle_360()

    # Try common Execute signatures
    try:
        res = RevolveFaces.Execute(FaceSelection.Create(face), axis_line, ang, opts)
    except:
        # some builds swap order or take different selection wrappers
        res = RevolveFaces.Execute(Selection.Create(face), axis_line, ang, opts)

    created = list(res.CreatedBodies)
    if not created:
        raise Exception("RevolveFaces returned no bodies (profile may be open or self-intersecting).")

    tool = created[0]
    tool.SetName(name)

    # Delete the temp sketch-fill body (optional cleanup)
    try:
        Delete.Execute(Selection.Create(temp))
    except:
        pass

    return tool

def create_two_pulley_rig():
    comp = get_or_create_component("RigidParts_2Pulley")

    W_ow = W_outer + 2.0*t_overwrap
    H_ow = H_outer + 2.0*t_overwrap

    zc = 0.5 * L_extrude

    # Center of pulleys above/below the cable mid-height
    y_offset = 0.5*H_ow + Pulley_Gap + Pulley_R

    # -------------------------
    # Helper: create wheel body
    # -------------------------
    def create_wheel(name, y0):
        set_sketch_plane_yz_at_x(-0.5 * Pulley_Width)
        sketch_circle(y0, zc, Pulley_R)
        wheel = extrude_last_profile_along_x(name, Pulley_Width, pick_largest=True)
        return wheel

    # --------------------------------------------
    # Helper: build + place + subtract groove tool
    # --------------------------------------------
    def cut_groove_into_wheel(wheel, tool_name, y0):
        # Build tool at x=0, centered at (Y=pulley_R, Z=0) by construction
        tool = create_groove_tool_revolved_about_x(
            tool_name,
            pulley_R=Pulley_R,
            W_ow=W_ow,
            H_ow=H_ow,
            clear=Pulley_Clear
        )

        # Move tool so its center matches wheel center:
        # tool currently centered at (Y=Pulley_R, Z=0). We need (Y=y0, Z=zc).
        dy = y0 - Pulley_R
        dz = zc - 0.0

        Move.Translate(Selection.Create(tool), Direction.DirY, MM(dy))
        Move.Translate(Selection.Create(tool), Direction.DirZ, MM(dz))

        # Subtract
        subtract_tool_from_target(wheel, tool)

        # Cleanup tool body
        try:
            Delete.Execute(Selection.Create(tool))
        except:
            pass

    # ==============
    # TOP pulley
    # ==============
    wheel_top = create_wheel("Rig_Pulley_Top", y0=+y_offset)
    cut_groove_into_wheel(wheel_top, "TMP_GROOVE_Top", y0=+y_offset)
    wheel_top = move_body_to_component(wheel_top, comp)
    create_ns("Rig_Pulley_Top", wheel_top)

    # ==============
    # BOTTOM pulley
    # ==============
    wheel_bot = create_wheel("Rig_Pulley_Bot", y0=-y_offset)
    cut_groove_into_wheel(wheel_bot, "TMP_GROOVE_Bot", y0=-y_offset)
    wheel_bot = move_body_to_component(wheel_bot, comp)
    create_ns("Rig_Pulley_Bot", wheel_bot)

    return wheel_top, wheel_bot
#==========================================================


# -----------------------------
# 5. FACE-LEVEL NAMED SELECTIONS
# -----------------------------
def get_body(name):
    for b in GetRootPart().Bodies:
        if b.Name == name:
            return b
    raise Exception("Body not found: " + name)
#shield = find_body_anywhere_by_name("Shield")[0]
#overwrap = find_body_anywhere_by_name("Overwrap")[0]
#second_extrusion = find_body_anywhere_by_name("Second_Extrusion")[0]
#create_layer_side_ns_by_normal(shield, "NS_Shield_Outer", "NS_Shield_Inner")
#create_layer_side_ns_by_normal(overwrap, "NS_Overwrap_Outer", "NS_Overwrap_Inner")
#create_layer_side_ns_by_normal(second_extrusion, "NS_Second_Extrusion_Outer", "NS_Second_Extrusion_Inner")


################# Stiffness Calculations ############################


def get_body_stiffness_contribution(body, E_modulus):
    """
    Calculates EI contribution using Parallel Axis Theorem.
    Uses Face properties to ensure 2D cross-section accuracy.
    """
    if body is None or not is_valid(body):
        return 0.0, 0.0

    # Locate the cross-section face (XY plane)
    target_face = _largest_xy_face(body)
    if target_face is None:
        return 0.0, 0.0

    # Get MassProperties of the FACE
    # In SC API, for a face, 'Mass' often returns the Area
    face_sel = FaceSelection.Create([target_face])
    mp = MeasureHelper.GetMassProperties(face_sel)
    
    # Robustly handle Area and Centroid from the MP object
    # If .Area isn't found, .Mass is the fallback for 2D objects
    try:
        area = float(mp.Area)
    except:
        area = float(mp.Mass) 
        
    # Area Centroid (d values for Parallel Axis Theorem)
    # Distance to Global Y-axis (x displacement) and Global X-axis (y displacement)
    cx = mp.Centroid.X
    cy = mp.Centroid.Y

    # Second Moment of Inertia about the face's own centroid
    # In V252, mp.Inertia returns the Moments of Inertia tensor
    I_local_x = mp.Inertia.XX
    I_local_y = mp.Inertia.YY

    # Parallel Axis Theorem: I_global = I_local + (Area * distance^2)
    # EI_x uses Y-offset (cy); EI_y uses X-offset (cx)
    I_global_x = I_local_x + (area * (cy**2))
    I_global_y = I_local_y + (area * (cx**2))

    return E_modulus * I_global_x, E_modulus * I_global_y

def pick_smallest_face(body):
    faces = []
    for f in list(body.Faces):
        try:
            faces.append((float(f.Area), f))
        except:
            pass
    faces.sort(key=lambda t: t[0])
    return faces[0][1] if faces else None


def _print_tensor(label, tensor, scale=1e12):
    # tensor assumed in m^4; scale to mm^4 by default
    try:
        m = [[tensor.GetValue(i, j) * scale for j in range(3)] for i in range(3)]
        print(label)
        print("  [{: .6e} {: .6e} {: .6e}]".format(m[0][0], m[0][1], m[0][2]))
        print("  [{: .6e} {: .6e} {: .6e}]".format(m[1][0], m[1][1], m[1][2]))
        print("  [{: .6e} {: .6e} {: .6e}]".format(m[2][0], m[2][1], m[2][2]))
    except Exception as e:
        print(label, "  (failed to print tensor)", e)

def debug_inertia_for_body_face(body, face):
    print("\n" + "="*110)
    print("DEBUG BODY:", getattr(body, "Name", "<unnamed>"))
    print("Face area (mm^2):", float(face.Area) * 1e6)

    # Two selection pathways for the *same face*
    sel_generic = None
    try:
        sel_generic = Selection.Create(face)
    except Exception as e:
        print("Selection.Create(face) failed:", e)

    sel_face = None
    try:
        sel_face = FaceSelection.Create([face])
    except Exception as e:
        print("FaceSelection.Create([face]) failed:", e)

    for sel_name, sel in [("Selection.Create(face)", sel_generic),
                          ("FaceSelection.Create([face])", sel_face)]:
        if sel is None:
            continue

        print("\n---", sel_name, "---")
        try:
            c = MeasureHelper.GetCentroid(sel)
            print("Centroid (mm):", (float(c.X)*1000.0, float(c.Y)*1000.0, float(c.Z)*1000.0))
        except Exception as e:
            print("GetCentroid failed:", e)
            c = None

        try:
            mp = MeasureHelper.GetMassProperties(sel)
        except Exception as e:
            print("GetMassProperties failed:", e)
            continue

        # Tell us whether SC thinks this is a SOLID vs a LAMINA-ish thing
        for attr in ["Area", "Volume", "Mass"]:
            try:
                print("mp.{:>6s} =".format(attr), getattr(mp, attr))
            except:
                pass

        # Frames: origin + centroid
        origin = Point.Create(0, 0, 0)
        frame_O = Frame.Create(origin, Direction.DirX, Direction.DirY)

        if c is None:
            # if centroid wasn't available, still print origin-based tensor
            c = origin
        frame_C = Frame.Create(c, Direction.DirX, Direction.DirY)

        # Print tensors about origin and about centroid
        try:
            T_O = mp.GetInertiaTensor(frame_O)
            _print_tensor("Inertia tensor about ORIGIN (X,Y,Z) [mm^4]:", T_O, scale=1e12)
        except Exception as e:
            print("GetInertiaTensor(origin frame) failed:", e)

        try:
            T_C = mp.GetInertiaTensor(frame_C)
            _print_tensor("Inertia tensor about CENTROID (X,Y,Z) [mm^4]:", T_C, scale=1e12)
        except Exception as e:
            print("GetInertiaTensor(centroid frame) failed:", e)

        # Also print just diagonals clearly
        try:
            T_C = mp.GetInertiaTensor(frame_C)
            Ixx = T_C.GetValue(0,0) * 1e12
            Iyy = T_C.GetValue(1,1) * 1e12
            Izz = T_C.GetValue(2,2) * 1e12
            Ixy = T_C.GetValue(0,1) * 1e12
            print("Centroid diag [mm^4]: Ixx=", Ixx, " Iyy=", Iyy, " Izz=", Izz, "  Ixy=", Ixy)
        except:
            pass

    print("="*110 + "\n")

# --- EXECUTION ---

def poly_area_centroid_I(points_xy):
    """
    points_xy: list of (x,y) in *mm*, NOT necessarily explicitly closed.
    Returns:
      A (mm^2), Cx (mm), Cy (mm),
      Ixx_O (mm^4) about global X-axis through origin,
      Iyy_O (mm^4) about global Y-axis through origin
    Signed area convention (CCW positive).
    """
    pts = points_xy[:]
    if len(pts) < 3:
        return 0.0, 0.0, 0.0, 0.0, 0.0

    # ensure closed
    if (pts[0][0] != pts[-1][0]) or (pts[0][1] != pts[-1][1]):
        pts.append(pts[0])

    A2 = 0.0  # 2A
    Cx6A = 0.0
    Cy6A = 0.0
    Ixx12 = 0.0
    Iyy12 = 0.0

    for i in range(len(pts) - 1):
        x0, y0 = pts[i]
        x1, y1 = pts[i+1]
        cross = x0*y1 - x1*y0
        A2 += cross
        Cx6A += (x0 + x1) * cross
        Cy6A += (y0 + y1) * cross
        Ixx12 += (y0*y0 + y0*y1 + y1*y1) * cross
        Iyy12 += (x0*x0 + x0*x1 + x1*x1) * cross

    A = 0.5 * A2
    if abs(A) < 1e-18:
        return 0.0, 0.0, 0.0, 0.0, 0.0

    Cx = Cx6A / (6.0 * A)
    Cy = Cy6A / (6.0 * A)
    Ixx_O = Ixx12 / 12.0
    Iyy_O = Iyy12 / 12.0
    return A, Cx, Cy, Ixx_O, Iyy_O


def _pt_key(p, tol=1e-9):
    # key in model units (meters); stable enough for connectivity
    return (round(float(p.X)/tol)*tol, round(float(p.Y)/tol)*tol, round(float(p.Z)/tol)*tol)

def _get_curve_obj_from_edge(edge):
    """
    Different SpaceClaim versions hang the curve/evaluator off different attributes.
    We probe a few common ones.
    """
    for attr in ("Geometry", "Curve", "Shape", "Underlying", "Definition"):
        try:
            obj = getattr(edge, attr)
            if obj is not None:
                return obj
        except:
            pass
    return None

def _curve_endpoints(curve_obj):
    """
    Get endpoints from curve object without assuming Evaluator exists.
    """
    # 1) If it has Evaluator
    try:
        ev = curve_obj.Evaluator
        tmin, tmax = ev.ParamRange
        p0 = ev.Evaluate(tmin).Point
        p1 = ev.Evaluate(tmax).Point
        return p0, p1, ev, (tmin, tmax)
    except:
        pass

    # 2) If it has StartPoint/EndPoint
    try:
        p0 = curve_obj.StartPoint
        p1 = curve_obj.EndPoint
        return p0, p1, None, None
    except:
        pass

    # 3) If it has GetEndPoints()
    try:
        p0, p1 = curve_obj.GetEndPoints()
        return p0, p1, None, None
    except:
        pass

    raise Exception("Could not extract endpoints from curve object: {}".format(type(curve_obj)))

def _edge_sample_points(edge, forward=True, n=40):
    """
    Returns list of Point objects in model units.
    Tries to sample via curve evaluator when available; falls back to endpoints.
    """
    curve_obj = _get_curve_obj_from_edge(edge)
    if curve_obj is None:
        raise Exception("Edge has no accessible curve/geometry: {}".format(type(edge)))

    p0, p1, ev, tr = _curve_endpoints(curve_obj)

    if ev is None or tr is None:
        pts = [p0, p1]
        if not forward:
            pts.reverse()
        return pts

    tmin, tmax = tr
    ts = [tmin + (tmax - tmin) * (i / float(n)) for i in range(n + 1)]
    pts = [ev.Evaluate(t).Point for t in ts]

    if not forward:
        pts.reverse()
    return pts

def _edge_endpoints(edge):
    """
    Uses the sampler to get endpoints (robust).
    """
    pts = _edge_sample_points(edge, forward=True, n=1)
    return pts[0], pts[-1]
def _edge_sample_points(edge, forward=True, n=30):
    """
    Samples points along an edge in model units (meters), returns list of Point objects.
    forward=True means from v0->v1, else reversed.
    """
    try:
        ev = edge.Evaluator
        tmin, tmax = ev.ParamRange
        ts = [tmin + (tmax - tmin) * (i / float(n)) for i in range(n+1)]
        pts = [ev.Evaluate(t).Point for t in ts]
    except:
        # fallback: just endpoints
        p0, p1 = _edge_endpoints(edge)
        pts = [p0, p1]

    if not forward:
        pts.reverse()
    return pts

def _build_loops_from_edges(edges):
    """
    Returns loops as list of (edge, forward_bool) sequences.
    """
    # Map vertex -> incident edges
    from collections import defaultdict
    v2e = defaultdict(list)
    edge_info = []
    for e in edges:
        p0, p1 = _edge_endpoints(e)
        k0, k1 = _pt_key(p0), _pt_key(p1)
        edge_info.append((e, k0, k1))
        v2e[k0].append(e)
        v2e[k1].append(e)

    unused = set(edges)
    loops = []

    # helper to get the "other" vertex key for edge, given current vertex key
    def other_key(e, curk):
        for (ee, k0, k1) in edge_info:
            if ee is e:
                return k1 if curk == k0 else k0
        return None

    # helper: does edge go forward from curk?
    def edge_forward_from(e, curk):
        for (ee, k0, k1) in edge_info:
            if ee is e:
                return (curk == k0)
        return True

    while unused:
        e0 = next(iter(unused))
        unused.remove(e0)

        # start at one endpoint
        p0, p1 = _edge_endpoints(e0)
        startk = _pt_key(p0)
        curk = _pt_key(p1)

        loop = [(e0, True)]  # e0 forward p0->p1
        safety = 0
        while curk != startk and safety < 10000:
            safety += 1
            # pick next unused edge incident to curk
            candidates = [e for e in v2e[curk] if e in unused]
            if not candidates:
                break

            enext = candidates[0]
            unused.remove(enext)

            fwd = edge_forward_from(enext, curk)
            loop.append((enext, fwd))
            curk = other_key(enext, curk)

        loops.append(loop)

    return loops

def face_boundary_loops_as_polylines_xy_mm(face, samples_per_edge=30):
    """
    Returns list of loops, each loop is list of (x_mm, y_mm) points.
    """
    # collect boundary edges of the face
    edges = list(face.Edges)
    loops = _build_loops_from_edges(edges)

    out = []
    for loop in loops:
        pts_xy = []
        for (e, fwd) in loop:
            pts = _edge_sample_points(e, forward=fwd, n=samples_per_edge)
            # append, avoiding duplicate join point
            for j, p in enumerate(pts):
                x_mm = float(p.X) * 1000.0
                y_mm = float(p.Y) * 1000.0
                if pts_xy and j == 0:
                    # skip first point if it matches last
                    if abs(x_mm - pts_xy[-1][0]) < 1e-6 and abs(y_mm - pts_xy[-1][1]) < 1e-6:
                        continue
                pts_xy.append((x_mm, y_mm))
        out.append(pts_xy)
    return out


def section_props_from_face_polygon(face, samples_per_edge=30):
    """
    Returns:
      A_mm2, Cx_mm, Cy_mm, Ixx_centroid_mm4, Iyy_centroid_mm4
    """
    loops = face_boundary_loops_as_polylines_xy_mm(face, samples_per_edge=samples_per_edge)
    if not loops:
        return 0.0, 0.0, 0.0, 0.0, 0.0

    # compute per-loop properties about origin
    per = []
    for pts in loops:
        A, Cx, Cy, Ixx_O, Iyy_O = poly_area_centroid_I(pts)
        per.append((A, Cx, Cy, Ixx_O, Iyy_O, pts))

    # pick outer as largest magnitude area
    per.sort(key=lambda t: abs(t[0]), reverse=True)
    outer = per[0]
    holes = per[1:]

    # accumulate origin-based totals with sign: outer + holes (holes should subtract)
    # If a hole comes in CCW positive, force it negative (and vice versa).
    def force_negative(A, Cx, Cy, Ixx_O, Iyy_O):
        if A > 0:
            return (-A, Cx, Cy, -Ixx_O, -Iyy_O)
        return (A, Cx, Cy, Ixx_O, Iyy_O)

    A_tot, Cx_num, Cy_num, Ixx_O_tot, Iyy_O_tot = 0.0, 0.0, 0.0, 0.0, 0.0

    # outer: keep its sign as-is
    Ao, Cxo, Cyo, Ixxo, Iyyo, _ = outer
    A_tot += Ao
    Cx_num += Cxo * Ao
    Cy_num += Cyo * Ao
    Ixx_O_tot += Ixxo
    Iyy_O_tot += Iyyo

    # holes: force negative contribution
    for (A, Cx, Cy, Ixx_O, Iyy_O, _) in holes:
        A2, Cx2, Cy2, Ixx2, Iyy2 = force_negative(A, Cx, Cy, Ixx_O, Iyy_O)
        A_tot += A2
        Cx_num += Cx2 * A2
        Cy_num += Cy2 * A2
        Ixx_O_tot += Ixx2
        Iyy_O_tot += Iyy2

    if abs(A_tot) < 1e-18:
        return 0.0, 0.0, 0.0, 0.0, 0.0

    Cx_tot = Cx_num / A_tot
    Cy_tot = Cy_num / A_tot

    # shift to centroidal moments
    Ixx_C = Ixx_O_tot - A_tot * (Cy_tot**2)
    Iyy_C = Iyy_O_tot - A_tot * (Cx_tot**2)

    return abs(A_tot), Cx_tot, Cy_tot, Ixx_C, Iyy_C

def calculate_explicit_composite_stiffness(modulus_map):
    target_names = ["conductor[1]", "conductor[2]", "Shield", "Overwrap", "Second_Extrusion"]
    if core_mode == "separate":
        target_names += ["single_core[1]", "single_core[2]"]
    else:
        target_names += ["single_core_merged"]
    if drain_opt:
        target_names.extend(["drain[1]", "drain[2]"])
        
    total_EIx = 0.0
    total_EIy = 0.0
    
    print("\n" + "="*95)
    print("{:<20} | {:<10} | {:<10} | {:<12} | {:<12}".format(
        "Body Name", "Area(mm2)", "Cx(mm)", "Ix ", "Iy "))
    print("-" * 95)

    for name in target_names:
        body, _ = find_body_anywhere_by_name(name)
        if body is None: continue

        E_val = 0.0
        for key, value in modulus_map.items():
            if key in name:
                E_val = value
                break
        if E_val <= 0: continue

        # 1. Pick the correct XY face
        faces_list = sorted([(float(f.Area), f) for f in body.Faces], key=lambda t: t[0])
        if not faces_list: continue
        face = faces_list[0][1]

        # 2. FORCE FaceSelection - DO NOT use Selection.Create(face)
        # This is the critical step to avoid the 3D-slab error.
        f_sel = FaceSelection.Create([face])
        
      # 1. Get the physical length of the body in METERS
        # This is the "thickness" that is inflating your mass properties
       # L_m = get_body_length_z_m(body) 
        
        # Use polygon/Green's theorem on the face boundary
        A_mm2, cx_mm, cy_mm, Ixx_c_mm4, Iyy_c_mm4 = section_props_from_face_polygon(face, samples_per_edge=40)

        # If you want EI about global origin (0,0): parallel axis from centroidal to origin
        Ixx_O_mm4 = Ixx_c_mm4 + A_mm2 * (cy_mm**2)
        Iyy_O_mm4 = Iyy_c_mm4 + A_mm2 * (cx_mm**2)

        body_eix = E_val * Ixx_O_mm4
        body_eiy = E_val * Iyy_O_mm4
                
        total_EIx += body_eix
        total_EIy += body_eiy
        
        print("{:<20} | {:<10.4f} | {:<10.4f} | {:<12.4e} | {:<12.4e}".format(
            name, area_m2*1e6, c_pt.X*1000.0, body_eix, body_eiy))

    print("="*95)
    print("TOTAL COMPOSITE EIx: {:.6e} N-mm2".format(total_EIx))
    print("TOTAL COMPOSITE EIy: {:.6e} N-mm2".format(total_EIy))
    return total_EIx, total_EIy

def get_polygon_props(face):
    """
    Calculates Area, Ix, and Iy using Green's Theorem by sampling face loops.
    Accesses the underlying Modeler Face to find Loops.
    """
    total_area = 0.0
    total_ixx = 0.0
    total_iyy = 0.0

    # Access the underlying geometry shape to get Loops
    modeler_face = face.Shape 

    for loop in modeler_face.Loops:
        # Get a faceted approximation of the loop
        # Tighten MM(0.001) if you need higher precision
        tess = loop.GetTessellation(MM(0.001), 0)
        points = [(pt.X * 1000.0, pt.Y * 1000.0) for pt in tess]

        # Standard Polygon Integration (Green's Theorem)
        area = 0.0
        ixx = 0.0
        iyy = 0.0
        n = len(points)
        
        if n < 3: continue

        for i in range(n):
            x0, y0 = points[i]
            x1, y1 = points[(i + 1) % n]
            
            cross = x0 * y1 - x1 * y0
            area += cross
            ixx += (y0**2 + y0*y1 + y1**2) * cross
            iyy += (x0**2 + x0*x1 + x1**2) * cross
            
        area /= 2.0
        ixx /= 12.0
        iyy /= 12.0

        # SpaceClaim Loops: Outer is CCW (positive), Inner is CW (negative).
        # Summing them naturally subtracts holes from the total.
        total_area += area
        total_ixx += ixx
        total_iyy += iyy
    
    return abs(total_area), abs(total_ixx), abs(total_iyy)

def calculate_explicit_composite_stiffness_v2(modulus_map):
    target_names = ["conductor[1]", "conductor[2]", "Shield", "Overwrap", "Second_Extrusion"]
    if core_mode == "separate":
        target_names += ["single_core[1]", "single_core[2]"]
    else:
        target_names += ["single_core_merged"]
        
    if drain_opt:
        target_names.extend(["drain[1]", "drain[2]"])
        
    total_EIx = 0.0
    total_EIy = 0.0
    
    print("\n" + "="*105)
    print("{:<20} | {:<12} | {:<12} | {:<15} | {:<15}".format(
        "Body Name", "Area (mm2)", "Cx (mm)", "EIx (N-mm2)", "EIy (N-mm2)"))
    print("-" * 105)

    for name in target_names:
        body, _ = find_body_anywhere_by_name(name)
        if body is None:
            continue

        # Match E-modulus (MPa)
        E_val = 0.0
        for key, value in modulus_map.items():
            if key in name:
                E_val = value
                break
        
        if E_val <= 0:
            continue

        # Pick the largest face (the cross-section)
        # We use face.Area (meters) for sorting
        faces_list = sorted([(float(f.Area), f) for f in body.Faces], key=lambda t: t[0], reverse=True)
        if not faces_list:
            continue
        face = faces_list[0][1]

        # Calculate geometric properties via boundary sampling (Pure Geometry)
        area_mm2, Ix_mm4, Iy_mm4 = get_polygon_props(face)
        
        # Calculate stiffness contributions (MPa * mm^4 = N-mm^2)
        body_eix = E_val * Ix_mm4
        body_eiy = E_val * Iy_mm4
        
        total_EIx += body_eix
        total_EIy += body_eiy
        
        # Determine centroid for reporting (using World Centroid)
        f_sel = FaceSelection.Create([face])
        centroid_pt = MeasureHelper.GetCentroid(f_sel)
        
        print("{:<20} | {:<12.4f} | {:<12.4f} | {:<15.4e} | {:<15.4e}".format(
            name, area_mm2, centroid_pt.X * 1000.0, body_eix, body_eiy))

    print("="*105)
    print("TOTAL COMPOSITE EIx: {:.6e} N-mm2".format(total_EIx))
    print("TOTAL COMPOSITE EIy: {:.6e} N-mm2".format(total_EIy))
    
    return total_EIx, total_EIy

material_moduli = {
    "conductor": 110000.0,
    "drain": 110000.0,
    "Shield": 3000.0,
    "single_core": 500.0,
    "Second_Extrusion": 25.0,
    "Overwrap": 3650.0
}

if run_stiffness:
    calculate_explicit_composite_stiffness(material_moduli)
    
    ###################
    # --- units helper (keep if MM isn't already defined in your script) ---
try:
    MM(1.0)
except:
    def MM(x):  # mm -> meters (SC internal)
        return x / 1000.0

def P2(x, y):
    return Point2D.Create(MM(x), MM(y))

def solidify_sketch():
    ViewHelper.SetViewMode(InteractionMode.Solid)

def set_sketch_plane_yz_at_x(x0):
    frame = Frame.Create(
        Point.Create(MM(x0), MM(0), MM(0)),
        Direction.DirY,
        Direction.DirZ
    )
    ViewHelper.SetSketchPlane(Plane.Create(frame))

def sketch_circle(cu, cv, r):
    # In a YZ sketch, Point2D(x,y) corresponds to (Y,Z).
    return SketchCircle.Create(Point2D.Create(MM(cu), MM(cv)), MM(r)).CreatedCurves

def pick_planar_face(body):
    faces = list(body.Faces)
    if not faces:
        raise Exception("No faces found (sketch not closed?)")
    # largest planar face is usually correct
    return max(faces, key=lambda f: f.Area)

def angle_360():
    # robust across builds
    try:
        return Angle.Create(2.0 * math.pi)  # radians
    except:
        try:
            return Angle.Degrees(360.0)
        except:
            return 360.0

def extrude_last_profile_along_x(name, length_mm, pick_largest=True):
    solidify_sketch()
    temp_body = GetRootPart().Bodies[-1]
    face = max(list(temp_body.Faces), key=lambda f: f.Area) if pick_largest else list(temp_body.Faces)[0]
    opts = ExtrudeFaceOptions()
    opts.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirX, MM(length_mm), opts)
    created = list(res.CreatedBodies)
    if not created:
        raise Exception("Extrude created no bodies for '%s'" % name)
    b = created[0]
    b.SetName(name)
    return b

def subtract_tool_from_target(target_body, tool_body):
    # Subtract tool from target. Keep target; discard tool later.
    try:
        # Some builds accept (targets, tools)
        res = Boolean.Subtract(Selection.Create([target_body]), Selection.Create([tool_body]))
        return res
    except:
        # Fallback signature used in some versions
        res = Boolean.Subtract(Selection.Create([target_body, tool_body]))
        return res

# --------- IMPORTANT: You need this function if you use Double-D groove ---------
def sketch_doubleD_shifted_YZ(y_center, z_center, dx, R):
    # Draw CLOSED Double-D loop in the current YZ sketch.
    SketchArc.CreateSweepArc(P2(y_center + dx, z_center),
                             P2(y_center + dx, z_center + R),
                             P2(y_center + dx, z_center - R),
                             True)
    SketchArc.CreateSweepArc(P2(y_center - dx, z_center),
                             P2(y_center - dx, z_center + R),
                             P2(y_center - dx, z_center - R),
                             False)
    SketchLine.Create(P2(y_center - dx, z_center + R), P2(y_center + dx, z_center + R))
    SketchLine.Create(P2(y_center - dx, z_center - R), P2(y_center + dx, z_center - R))

def create_groove_tool_revolved_about_x(name, pulley_R, W_ow, H_ow, clear):
    # Groove tool is a torus-like body created by revolving a Double-D profile about X axis.
    Wg = W_ow + 2.0*clear
    Hg = H_ow + 2.0*clear
    Rg  = 0.5*Hg
    dxg = 0.5*(Wg - Hg)

    yc = pulley_R
    zc = 0.0

    set_sketch_plane_yz_at_x(0.0)
    sketch_doubleD_shifted_YZ(yc, zc, dxg, Rg)

    solidify_sketch()
    temp = GetRootPart().Bodies[-1]
    face = pick_planar_face(temp)

    axis_line = Line.Create(Point.Create(MM(0), MM(0), MM(0)), Direction.DirX)
    opts = RevolveFaceOptions()
    opts.RevolveType = RevolveType.ForceIndependent

    ang = angle_360()
    try:
        res = RevolveFaces.Execute(FaceSelection.Create(face), axis_line, ang, opts)
    except:
        res = RevolveFaces.Execute(Selection.Create(face), axis_line, ang, opts)

    created = list(res.CreatedBodies)
    if not created:
        raise Exception("RevolveFaces returned no bodies (profile may be open/self-intersecting).")

    tool = created[0]
    tool.SetName(name)

    # cleanup temp sketch-fill
    try:
        Delete.Execute(Selection.Create(temp))
    except:
        pass

    return tool

# ============================================================
# Two-pulley rig builder
# ============================================================

# You can expose these as Parameters[...] if you want; here are safe defaults:
Pulley_R      = 10.0   # mm (wheel radius)
Pulley_Width  = 6.0    # mm (wheel thickness along X)
Pulley_Gap    = 0.20   # mm (gap between cable OD and groove centerline placement)
Pulley_Clear  = 0.05   # mm (extra groove clearance)

def create_two_pulley_rig():
    W_ow = W_outer + 2.0*t_overwrap
    H_ow = H_outer + 2.0*t_overwrap

    zc = 0.5 * L_extrude
    y_offset = 0.5*H_ow + Pulley_Gap + Pulley_R

    # --- wheel (simple cylinder) ---
    def create_wheel(name, y0):
        set_sketch_plane_yz_at_x(-0.5 * Pulley_Width)
        sketch_circle(y0, zc, Pulley_R)  # in YZ plane: (Y,Z)
        wheel = extrude_last_profile_along_x(name, Pulley_Width, pick_largest=True)
        return wheel

    # --- cut groove ---
    def cut_groove_into_wheel(wheel, tool_name, y0):
        tool = create_groove_tool_revolved_about_x(
            tool_name,
            pulley_R=Pulley_R,
            W_ow=W_ow,
            H_ow=H_ow,
            clear=Pulley_Clear
        )

        # tool as-built is centered at (Y=Pulley_R, Z=0); move to (Y=y0, Z=zc)
        dy = y0 - Pulley_R
        dz = zc - 0.0
        Move.Translate(Selection.Create(tool), Direction.DirY, MM(dy))
        Move.Translate(Selection.Create(tool), Direction.DirZ, MM(dz))

        subtract_tool_from_target(wheel, tool)

        # delete tool
        try:
            Delete.Execute(Selection.Create(tool))
        except:
            pass

    # TOP pulley
    wheel_top = create_wheel("Rig_Pulley_Top", y0=+y_offset)
    cut_groove_into_wheel(wheel_top, "TMP_GROOVE_Top", y0=+y_offset)

    # BOTTOM pulley
    wheel_bot = create_wheel("Rig_Pulley_Bot", y0=-y_offset)
    cut_groove_into_wheel(wheel_bot, "TMP_GROOVE_Bot", y0=-y_offset)

    return wheel_top, wheel_bot

# ---- call it ----
# wheel_top, wheel_bot = create_two_pulley_rig()

# RUN PIPELINE
# ============================================================
if run_benchmark_rig:
    mode = bending_mode.lower().strip()
    if mode not in ["good_3point", "bad_3point", "good_2pulleys", "bad_2pulleys"]:
        raise Exception("bending_mode must be 'good_3point', 'bad_3point', 'good_2pulleys', or 'bad_2pulleys'")

    # 3-point rig
    if mode in ["good_3point", "bad_3point"]:
        nose_body = create_loading_nose(mode)
        s1, s2 = create_two_supports(mode)

        create_ns("Rig_Support_1", s1)
        create_ns("Rig_Support_2", s2)

        # optional cleanup of stray groups
        for k in range(1, 10):
            try:
                NamedSelection.Delete("Group%d" % k)
            except:
                pass

    # 2-pulley rig
    elif mode in ["good_2pulleys", "bad_2pulleys"]:
        wheel_top, wheel_bot = create_two_pulley_rig()

        # If your create_two_pulley_rig already does create_ns(), you can omit these.
        # Otherwise:
        try: create_ns("Rig_Pulley_Top", wheel_top)
        except: pass
        try: create_ns("Rig_Pulley_Bot", wheel_bot)
        except: pass

        for k in range(1, 10):
            try:
                NamedSelection.Delete("Group%d" % k)
            except:
                pass
        # Optional: named selections for supports if you want
    # create_ns(_body_name(s1), s1)
    # create_ns(_body_name(s2), s2)
