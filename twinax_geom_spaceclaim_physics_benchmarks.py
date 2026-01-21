# Python Script, API Version = V252


import math

from SpaceClaim.Api.V252 import *

from SpaceClaim.Api.V252.Geometry import *

from SpaceClaim.Api.V252.Modeler import *

ClearAll()


# -----------------------------

# 1. PARAMETERS / INPUTS

# -----------------------------

def get_param(name, default_val):

    try: return float(Parameters[name])

    except: return default_val


run_benchmark_rig = False  # Set to True to enable 3-point rig creation

# Logic options using the Parameters check pattern

filler_opt   = int(Parameters.filler_option) if hasattr(Parameters, 'filler_option') else 0
drain_opt    = float(Parameters.has_drains) if hasattr(Parameters, 'has_drains') else 0.0
is_a_doublet = int(Parameters.is_a_doublet) if hasattr(Parameters, 'is_a_doublet') else 0
is_elliptic  = int(Parameters.is_elliptic) if hasattr(Parameters, 'is_elliptic') else 1
h_mix        = float(Parameters.mixing_factor) if hasattr(Parameters, 'mixing_factor') else 0.7
bending_benchmark_opt = int(Parameters.bending_mode_option) if hasattr(Parameters, 'bending_mode_option') else 0
# filler_mode logic

filler_mode = "shell" if filler_opt == 1 else "fill"

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
C2C    = float(Parameters.c2c) if hasattr(Parameters, 'c2c') else 3.0
D_core = float(Parameters.diam_core) if hasattr(Parameters, 'diam_core') else 2.8

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


# Shield inner dimensions

W_in = W_outer - 2.0 * t_shield
H_in = H_outer - 2.0 * t_shield
R_in = H_in / 2.0
dx_in = (W_in - H_in) / 2.0
R_shield = H_outer / 2.0
dx_shield = (W_outer - H_outer) / 2.0


## Take care of the tangency situation between the two cores
overlap_tol = 0.012  # mm

# Geometry relationship between conductor spacing and core size

if (filler_mode == "shell"): 
    D_core = C2C
elif (filler_mode == "fill"):
    delta = C2C - D_core
    # snap away from zero-thickness tangency
    if abs(delta) < overlap_tol:
        # choose one:
        # C2C = D_core + overlap_tol            # Policy A (always spaced)
        C2C = D_core + overlap_tol if delta > 0 else D_core - overlap_tol  # Policy B
        delta = C2C - D_core

    core_mode = "spaced" if delta > 0 else "merged" 
dx_cond = C2C / 2.0
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

def share_topology_pair(body_a, body_b, verbose=True):
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

def extrude_and_name(name, length_mm, pick_largest=False):
    root = GetRootPart()
    before = list(root.Bodies)

    solidify_sketch()

    after = list(root.Bodies)
    new_bodies = [b for b in after if b not in before]

    if not new_bodies:
        raise Exception(
            "Solidify produced NO sketch-fill body for '%s'. "
            "Profile is likely not closed (tiny gap/overlap), or multiple regions." % name
        )

    temp_body = new_bodies[-1]

    face = _largest_xy_face(temp_body)
    if face is None:
        face = max(list(temp_body.Faces), key=lambda f: f.Area)

    options = ExtrudeFaceOptions()
    options.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirZ, MM(length_mm), options)

    created = list(res.CreatedBodies)
    if not created:
        raise Exception("Extrude produced no result for '%s' (wrong face or invalid profile)." % name)

    new_body = created[0]
    new_body.SetName(name)
    return new_body


def extrude_and_name_old(name, length_mm, pick_largest=False):
    solidify_sketch()
    temp_body = GetRootPart().Bodies[-1]

    face = _largest_xy_face(temp_body)
    if face is None:
        # fallback
        face = max(list(temp_body.Faces), key=lambda f: f.Area)

    options = ExtrudeFaceOptions()
    options.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirZ, MM(length_mm), options)

    new_body = list(res.CreatedBodies)[0]
    new_body.SetName(name)
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
# 3. BUILD GEOMETRY

# -----------------------------
# (1) Conductors

set_sketch_plane_xy()

sketch_circle(-dx_cond, 0, r_cond)

c1 = extrude_and_name("conductor[1]", L_extrude)

#create_ns("conductor[1]", c1)

set_sketch_plane_xy()

sketch_circle(dx_cond, 0, r_cond)

c2 = extrude_and_name("conductor[2]", L_extrude)

#create_ns("conductor[2]", c2)


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
        if core_mode == "spaced":
            sketch_core_with_lumens(-dx_cond)
            score1 = extrude_and_name("single_core[1]", L_extrude, True)

            sketch_core_with_lumens(dx_cond)
            score2 = extrude_and_name("single_core[2]", L_extrude, True)


        elif core_mode == "merged":
            # Build BOTH multilumen cores first (left + right), then Boolean-union them
            set_sketch_plane_xy()
            sketch_core_with_lumens(-dx_cond)
            core1 = extrude_and_name("temp_core1", L_extrude, True)

            sketch_core_with_lumens(+dx_cond)
            core2 = extrude_and_name("temp_core2", L_extrude, True)
            #score = extrude_and_name("single_core_merged", L_extrude, True)
            # Union
            score = union_bodies("single_core_merged", [core1, core2])
            #score = merge_bodies_keep_name("single_core_merged", [core1, core2])

            # only delete bodies that are NOT the output
            victims = []
            for b in [core1, core2]:
                if b is not score:
                    victims.append(b)
            if victims:
                try: Delete.Execute(BodySelection.Create(victims))
                except: pass

    else:

        if core_mode == "spaced":
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
            set_sketch_plane_xy()
            sketch_circle(-dx_cond, 0, r_core)
            sketch_circle(-dx_cond, 0, r_cond)
            core1 = extrude_and_name("temp_core1", L_extrude, True)

            set_sketch_plane_xy()
            sketch_circle(dx_cond, 0, r_core)
            sketch_circle(dx_cond, 0, r_cond)
            core2 = extrude_and_name("temp_core2", L_extrude, True)

            # Perform Boolean Merge (Union)
            score = union_bodies("single_core_merged", [core1, core2])

            # Assign final name and NS
            #merged_result.SetName("single_core_merged")
            #create_ns("single_core_merged", merged_result)
            

# (3) Filler / Doublet

set_sketch_plane_xy()

if is_a_doublet:

    sketch_profile(dx_in, R_in, W_in, H_in)
    sketch_circle(-dx_cond, 0, r_cond)
    sketch_circle(dx_cond, 0, r_cond)
    second_extrusion = extrude_and_name("Second_Extrusion", L_extrude, True)
    #create_ns("Second_Extrusion", second_extrusion)

elif filler_mode == "fill":

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

    sketch_circle(-dx_cond, 0, r_core)
    sketch_circle(dx_cond, 0, r_core)
    second_extrusion = extrude_and_name("Second_Extrusion", L_extrude, True)
    #create_ns("Second_Extrusion", second_extrusion)


# --- SHARE TOPOLOGY / IMPRINT CHAIN (do this BEFORE Named Selections) ---
# conductor ↔ core
#share_topology_pair(c1, score1)
#share_topology_pair(c2, score2)

# core ↔ second extrusion
#share_topology_pair(score1, second_extrusion)
#share_topology_pair(score2, second_extrusion)

# (optional) if you actually need shield interfaces conformal:
# share_topology_pair(shield, second_extrusion)
# share_topology_pair(shield, overwrap)


# (4) Shield
print("Shield params: W_outer=%.4f H_outer=%.4f dx_shield=%.4f R_shield=%.4f"
      % (W_outer, H_outer, dx_shield, R_shield))
set_sketch_plane_xy()
sketch_profile(dx_shield, R_shield, W_outer, H_outer)
shield = extrude_and_name("Shield", L_extrude, True)
#create_ns("Shield", shield)


# (5) Drains

x_drain = dx_shield + R_shield + r_drain

if drain_opt:

    set_sketch_plane_xy()
    sketch_circle(-x_drain, 0, r_drain)
    drain1 = extrude_and_name("drain[1]", L_extrude)
    #create_ns("drain[1]", drain1)
    set_sketch_plane_xy()
    sketch_circle(x_drain, 0, r_drain)
    drain2 = extrude_and_name("drain[2]", L_extrude)
    #create_ns("drain[2]", drain2)


# (6) Overwrap

if drain_opt:

    # Existing drain-wrap logic relies on Double-D tangency; 

    # This remains as previously defined

    set_sketch_plane_xy()
    sketch_wrap_loop_with_drains(0)
    sketch_wrap_loop_with_drains(t_overwrap)
    overwrap = extrude_and_name("Overwrap", L_extrude, True)
    #create_ns("Overwrap", overwrap)

else:

    # 1. Create the Outer Solid for Overwrap
    set_sketch_plane_xy()
    sketch_profile(dx_shield, R_shield + t_overwrap, W_outer + 2*t_overwrap, H_outer + 2*t_overwrap)
    overwrap = extrude_and_name("Overwrap", L_extrude, True)
    #create_ns("Overwrap", overwrap)



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


def body_volume_mm3(body):
    """
    Returns volume in mm^3 if possible. If the API returns internal units,
    the value is still consistent for 'near zero' detection.
    """

    # --- Attempt 1: MassProperties on body (some versions expose this directly)
    try:
        mp = body.MassProperties
        v = float(mp.Volume)
        return v
    except:
        pass

    # --- Attempt 2: GetVolume / Volume property on body
    for attr in ["GetVolume", "Volume"]:
        try:
            v = getattr(body, attr)
            if callable(v):
                v = v()
            v = float(v)
            return v
        except:
            pass

    # --- Attempt 3: Mass properties via a Selection (common pattern)
    # Different builds expose different helpers; try a couple names.
    sel = BodySelection.Create([body])

    # 3a) MeasureHelper.GetMassProperties(...)
    try:
        mp = MeasureHelper.GetMassProperties(sel)
        v = float(mp.Volume)
        return v
    except:
        pass

    # 3b) MassProperties.Create(...) or similar
    try:
        mp = MassProperties.Create(sel)
        v = float(mp.Volume)
        return v
    except:
        pass

    # If all fail, return None (we'll report it)
    return None


def delete_zero_volume_bodies_in_component(comp_name, vol_tol=1e-6, verbose=True):
    """
    Deletes bodies in the given component whose volume <= vol_tol.
    Also prints a volume report for later use.
    """

    comp = _comp_by_name(comp_name)
    if comp is None:
        print("Component not found:", comp_name)
        return []

    bodies = list(comp.Content.Bodies)
    if not bodies:
        print("No bodies found in component:", comp_name)
        return []

    report = []
    victims = []
    unknown = []

    for b in bodies:
        v = body_volume_mm3(b)
        nm = getattr(b, "Name", "<?>")

        if v is None:
            unknown.append(b)
            if verbose:
                print("VOL=??   ", nm, "(could not compute)")
            continue

        report.append((nm, v))
        if verbose:
            print("VOL=%g  %s" % (v, nm))

        if abs(v) <= vol_tol:
            victims.append(b)

    if victims:
        Delete.Execute(BodySelection.Create(victims))
        print("Deleted %d zero-volume bodies (tol=%g) from %s" % (len(victims), vol_tol, comp_name))
    else:
        print("No zero-volume bodies found (tol=%g) in %s" % (vol_tol, comp_name))

    if unknown:
        print("WARNING: %d bodies had unknown volume (not deleted)." % len(unknown))

    # Return volume report sorted by absolute volume (useful later)
    report_sorted = sorted(report, key=lambda t: abs(t[1]))
    return report_sorted

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
if core_mode != "merged":
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
        cable_top_y = 0.5 * (H_outer + 2.0*t_overwrap)
        y_center = cable_top_y + Initial_Gap + r

        set_sketch_plane_yz_at_x(-0.5 * Nose_Length)
        sketch_circle(y_center, z_center, r)

        temp_nose = extrude_last_profile_along_x("Rig_Loading_Nose", Nose_Length)
        comp = get_or_create_component("RigidParts_3Point_Bending")
        final_nose = move_body_to_component(temp_nose, comp)
        create_ns("Rig_Loading_Nose", final_nose)
        return final_nose

    elif mode == "bad_3point":
        # Bad: nose on side (+X), sketch on ZX@y, extrude along Y
        cable_top_x = 0.5 * (W_outer + 2.0*t_overwrap + 2.0*D_drain if drain_opt else 0.0)
        x_center = cable_top_x + Initial_Gap + r

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
# RUN PIPELINE
# ============================================================


if run_benchmark_rig:
    mode = bending_mode.lower().strip()
    if mode not in ["good_3point", "bad_3point", "good_2pulleys", "bad_2pulleys"]:
        raise Exception("bending_mode must be 'good_3point', 'bad_3point', 'good_2pulleys', or 'bad_2pullyes'")



    # 3) Open supports (mode-dependent)
    if mode in ["good_3point", "bad_3point"]:
        nose_body = create_loading_nose(mode)
        s1, s2 = create_two_supports(mode)
        if mode == "good_3point":
           # s1 = split_by_ZX_faces_keep_bottom(s1, y_cut=0.0, final_name="Rig_Support_1")
            #s2 = split_by_ZX_faces_keep_bottom(s2, y_cut=0.0, final_name="Rig_Support_2")
            create_ns("Rig_Support_1", s1)
            create_ns("Rig_Support_2", s2)
        else:
           # s1 = split_by_YZ_faces_keep_bottom(s1, x_cut=0.0, final_name="Rig_Support_1")
           # s2 = split_by_YZ_faces_keep_bottom(s2, x_cut=0.0, final_name="Rig_Support_2")
            create_ns("Rig_Support_1", s1)
            create_ns("Rig_Support_2", s2)
        # Optional: named selections for supports if you want
    # create_ns(_body_name(s1), s1)
    # create_ns(_body_name(s2), s2)
    # Python Script, API Version = V252