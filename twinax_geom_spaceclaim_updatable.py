# Python Script, API Version = V22


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


# Logic options using the Parameters check pattern

filler_opt   = int(Parameters.filler_option) if hasattr(Parameters, 'filler_option') else 0
drain_opt    = float(Parameters.has_drains) if hasattr(Parameters, 'has_drains') else 0.0
is_a_doublet = int(Parameters.is_a_doublet) if hasattr(Parameters, 'is_a_doublet') else 0
is_elliptic  = int(Parameters.is_elliptic) if hasattr(Parameters, 'is_elliptic') else 1
h_mix        = float(Parameters.mixing_factor) if hasattr(Parameters, 'mixing_factor') else 0.7

# filler_mode logic

filler_mode = "shell" if filler_opt == 1 else "fill"

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
dx_cond = C2C / 2.0
if (filler_mode == "shell" and D_core != C2C): D_core = C2C

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



def create_ns(name, body_or_list):

    """Creates a Named Selection in the Groups tab for ANSYS Mechanical"""

    if isinstance(body_or_list, list):

        items = body_or_list

    else:

        items = [body_or_list]


    sel = BodySelection.Create(items)

    # Create the group and rename it

    design_group = NamedSelection.Create(sel, Selection.Empty())

    design_group.CreatedNamedSelection.SetName("NS_" + name)
  

# 3. BUILD GEOMETRY

# -----------------------------
# (1) Conductors

set_sketch_plane_xy()

sketch_circle(-dx_cond, 0, r_cond)

c1 = extrude_and_name("conductor[1]", L_extrude)

create_ns("conductor[1]", c1)

set_sketch_plane_xy()

sketch_circle(dx_cond, 0, r_cond)

c2 = extrude_and_name("conductor[2]", L_extrude)

create_ns("conductor[2]", c1)


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

        sketch_core_with_lumens(-dx_cond)
        score1 = extrude_and_name("single_core[1]", L_extrude, True)
        create_ns("single_core[1]", score1)

        # Build Core 2

        sketch_core_with_lumens(dx_cond)
        score2 = extrude_and_name("single_core[2]", L_extrude, True)
        create_ns("single_core[2]", score2)

    else:

        set_sketch_plane_xy()
        sketch_circle(-dx_cond, 0, r_core)
        sketch_circle(-dx_cond, 0, r_cond)
        score1 = extrude_and_name("single_core[1]", L_extrude, True)
        create_ns("single_core[1]", score1)

        set_sketch_plane_xy()
        sketch_circle(dx_cond, 0, r_core)
        sketch_circle(dx_cond, 0, r_cond)
        score2 = extrude_and_name("single_core[2]", L_extrude, True)
        create_ns("single_core[2]", score2)

# (3) Filler / Doublet

set_sketch_plane_xy()

if is_a_doublet:

    sketch_profile(dx_in, R_in, W_in, H_in)
    sketch_circle(-dx_cond, 0, r_cond)
    sketch_circle(dx_cond, 0, r_cond)
    second_extrusion = extrude_and_name("Second_Extrusion", L_extrude, True)
    create_ns("Second_Extrusion", second_extrusion)

elif filler_mode == "fill":

    sketch_profile(dx_in, R_in, W_in, H_in)
    sketch_circle(-dx_cond, 0, r_core)
    sketch_circle(dx_cond, 0, r_core)
    second_extrusion = extrude_and_name("Second_Extrusion", L_extrude, True)
    create_ns("Second_Extrusion", second_extrusion)

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
    create_ns("Second_Extrusion", second_extrusion)


# (4) Shield

set_sketch_plane_xy()
sketch_profile(dx_shield, R_shield, W_outer, H_outer)
shield = extrude_and_name("Shield", L_extrude, True)
create_ns("Shield", shield)


# (5) Drains

x_drain = dx_shield + R_shield + r_drain

if drain_opt:

    set_sketch_plane_xy()
    sketch_circle(-x_drain, 0, r_drain)
    drain1 = extrude_and_name("drain[1]", L_extrude)
    create_ns("drain[1]", drain1)
    set_sketch_plane_xy()
    sketch_circle(x_drain, 0, r_drain)
    drain2 = extrude_and_name("drain[2]", L_extrude)
    create_ns("drain[2]", drain2)


# (6) Overwrap

if drain_opt:

    # Existing drain-wrap logic relies on Double-D tangency; 

    # This remains as previously defined

    set_sketch_plane_xy()
    sketch_wrap_loop_with_drains(0)
    sketch_wrap_loop_with_drains(t_overwrap)
    overwrap = extrude_and_name("Overwrap", L_extrude, True)
    create_ns("Overwrap", overwrap)

else:

    # 1. Create the Outer Solid for Overwrap
    set_sketch_plane_xy()
    sketch_profile(dx_shield, R_shield + t_overwrap, W_outer + 2*t_overwrap, H_outer + 2*t_overwrap)
    overwrap = extrude_and_name("Overwrap", L_extrude, True)
    create_ns("Overwrap", overwrap)

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

# --- EXECUTE THE MOVE ---
move_all_root_bodies_to_component("cable_bodies")


# ------------------------------------------------------------
# ============================================================
# 7) RIGID 3-POINT BENDING: LOADING NOSE  +  TWO DOUBLE-D SUPPORTS
#   - Cable axis is Z (same as your extrusions)
#   - U/pocket profiles sketched in XY, extruded along Z
#   - Rigid bodies moved into component: RigidParts_3Point_Bending
#   - Supports located at 20% and 80% of L_extrude
# ============================================================

# -----------------------------
# Parameters (derived defaults)
# -----------------------------
Initial_Gap  = get_param("Initial_Gap", 0.05 * (H_outer + 2.0*t_overwrap))  # small gap above cable
Nose_Diam    = get_param("Nose_Diam", 2.0 * (H_outer + 2.0*t_overwrap))     # reasonable first pass
Nose_Length = get_param("Nose_Length", 4.0 * (H_outer + 2.0*t_overwrap))  # cylinder axis is X (extrude in X)

# -----------------------------
# Helpers (minimal + stable)
# -----------------------------
def set_sketch_plane_yz_at_x(x0):
    frame = Frame.Create(
        Point.Create(MM(x0), MM(0), MM(0)),
        Direction.DirY,
        Direction.DirZ
    )
    plane = Plane.Create(frame)
    ViewHelper.SetSketchPlane(plane)

def set_sketch_plane_xy_at_z(z0):
    frame = Frame.Create(
        Point.Create(MM(0), MM(0), MM(z0)),
        Direction.DirX,
        Direction.DirY
    )
    plane = Plane.Create(frame)
    ViewHelper.SetSketchPlane(plane)
    
def set_sketch_plane_zx_at_y(x0):
    frame = Frame.Create(
        Direction.DirX,
        Point.Create(MM(x0), MM(0), MM(0)),
        Direction.DirZ
    )
    plane = Plane.Create(frame)
    ViewHelper.SetSketchPlane(plane)

def get_or_create_component(name):
    root = GetRootPart()
    # Find existing
    for c in list(root.Components):
        try:
            if c.GetName() == name:
                return c
        except:
            if getattr(c, "Name", "") == name:
                return c

    # Create via Part blueprint + Component instance (your working style)
    doc = Window.ActiveWindow.Document
    part_def = Part.Create(doc, name)
    comp = Component.Create(root, part_def)
    return comp

def move_body_to_component_and_refetch(body, target_component, final_name):
    """
    Moves 'body' into component and returns a stable reference by refetching via temp name.
    This avoids 'Bodies[-1]' instability and prevents ping-pong / regen weirdness.
    """
    tmp_name = final_name + "__TEMP__"
    body.SetName(tmp_name)

    ComponentHelper.MoveBodiesToComponent(Selection.Create(body), target_component)

    for b in list(target_component.Content.Bodies):
        try:
            if b.GetName() == tmp_name:
                b.SetName(final_name)
                return b
        except:
            if getattr(b, "Name", "") == tmp_name:
                b.SetName(final_name)
                return b

    raise Exception("Moved body but could not refetch '%s' in target component." % final_name)

def extrude_largest_face_along_x(name, length_mm):
    solidify_sketch()
    temp_body = GetRootPart().Bodies[-1]
    faces = list(temp_body.Faces)
    if not faces:
        raise Exception("No faces found to extrude for %s (sketch not closed?)" % name)
    face = max(faces, key=lambda f: f.Area)

    options = ExtrudeFaceOptions()
    options.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirX, MM(length_mm), options)

    created = list(res.CreatedBodies)
    if not created:
        raise Exception("Extrude created no bodies for %s" % name)

    body = created[0]
    body.SetName(name)
    return body

def extrude_largest_face_along_z(name, length_mm):
    solidify_sketch()
    temp_body = GetRootPart().Bodies[-1]
    faces = list(temp_body.Faces)
    if not faces:
        raise Exception("No faces found to extrude for %s (sketch not closed?)" % name)
    face = max(faces, key=lambda f: f.Area)

    options = ExtrudeFaceOptions()
    options.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirZ, MM(length_mm), options)

    created = list(res.CreatedBodies)
    if not created:
        raise Exception("Extrude created no bodies for %s" % name)

    body = created[0]
    body.SetName(name)
    return body

def sketch_rectangle_xy(xmin, xmax, ymin, ymax):
    SketchLine.Create(P2(xmin, ymin), P2(xmax, ymin))
    SketchLine.Create(P2(xmax, ymin), P2(xmax, ymax))
    SketchLine.Create(P2(xmax, ymax), P2(xmin, ymax))
    SketchLine.Create(P2(xmin, ymax), P2(xmin, ymin))

def sketch_doubleD_shifted(x0, y0, dx, R):
    # right arc
    SketchArc.CreateSweepArc(P2(x0 + dx, y0), P2(x0 + dx, y0 + R), P2(x0 + dx, y0 - R), True)
    # left arc
    SketchArc.CreateSweepArc(P2(x0 - dx, y0), P2(x0 - dx, y0 + R), P2(x0 - dx, y0 - R), False)
    # bridges
    SketchLine.Create(P2(x0 - dx, y0 + R), P2(x0 + dx, y0 + R))
    SketchLine.Create(P2(x0 - dx, y0 - R), P2(x0 + dx, y0 - R))

# ============================================================
# 7A) Loading Nose (YZ circle -> X extrude, centered by construction)
def move_body_to_component(body, target_component):
    """
    Moves a body to a component and RETURNS THE NEW BODY.
    Crucial: The original 'body' reference is invalid after the move.
    """
    sel = Selection.Create(body)
    ComponentHelper.MoveBodiesToComponent(sel, target_component)
    
    # The body is now the last one in the target component's list
    # We must return this NEW reference.
    return target_component.Content.Bodies[-1]

def extrude_last_profile_along_x(name, length_mm, pick_largest=True, dir_sign=+1):
    solidify_sketch()
    temp_body = GetRootPart().Bodies[-1]
    face = max(temp_body.Faces, key=lambda f: f.Area) if pick_largest else temp_body.Faces[0]

    options = ExtrudeFaceOptions()
    options.ExtrudeType = ExtrudeType.ForceIndependent

    dirx = Direction.DirX if dir_sign > 0 else Direction.DirX.Negate()
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), dirx, MM(length_mm), options)
    new_body = res.CreatedBodies[0]
    new_body.SetName(name)
    return new_body

def create_loading_nose_in_rigid_component():
    r = Nose_Diam / 2.0

    # ... (Geometry Math is same) ...
    cable_top_y = 0.5 * (H_outer + 2.0 * t_overwrap)
    y_center = cable_top_y + Initial_Gap + r
    z_center = 0.5 * L_extrude

    # ... (Sketch & Extrude is same) ...
    set_sketch_plane_yz_at_x(-0.5 * Nose_Length)
    sketch_circle(y_center, z_center, r)
    
    # This 'nose' is temporary (lives in Root)
    temp_nose = extrude_last_profile_along_x("Rig_Loading_Nose", Nose_Length)

    # ... (Organize) ...
    rigid_comp = get_or_create_component("RigidParts_3Point_Bending")
    
    # FIX: Update 'nose' variable to the new body inside the component
    final_nose = move_body_to_component(temp_nose, rigid_comp)

    # ... (Named Selection) ...
    # Now we use 'final_nose', which is a valid reference
    create_ns("Rig_Loading_Nose", final_nose)

    return final_nose

# Execute
nose_body = create_loading_nose_in_rigid_component()

# ============================================================
# 7F) RIGID 3-POINT BENDING: TWO SUPPORTS WITH DOUBLE-D POCKET
#   - Pocket follows OUTER OVERWRAP shape (double-D envelope)
#   - Sketch in XY, extrude in Z (cable axis)
#   - Supports at 20% and 80% of L_extrude
#   - No boolean needed (pocket from sketch loops)
# ============================================================
def move_body_to_component_and_refetch(body, target_component, final_name):
    """
    Move body into component and refetch by name (stable).
    Avoids relying on 'Bodies[-1]' which can be unstable.
    """
    # Ensure unique temp name before move
    body.SetName(final_name + "__TEMP__")

    sel = Selection.Create(body)
    ComponentHelper.MoveBodiesToComponent(sel, target_component)

    # Refetch by name inside the target component
    for b in list(target_component.Content.Bodies):
        try:
            if b.GetName() == final_name + "__TEMP__":
                b.SetName(final_name)
                return b
        except:
            if getattr(b, "Name", "") == final_name + "__TEMP__":
                b.SetName(final_name)
                return b

    raise Exception("Moved body but could not refetch '%s' in target component." % final_name)

def set_sketch_plane_xy_at_z(z0):
    frame = Frame.Create(
        Point.Create(MM(0), MM(0), MM(z0)),
        Direction.DirX,
        Direction.DirY
    )
    plane = Plane.Create(frame)
    ViewHelper.SetSketchPlane(plane)

def sketch_rectangle_xy(xmin, xmax, ymin, ymax):
    SketchLine.Create(P2(xmin, ymin), P2(xmax, ymin))
    SketchLine.Create(P2(xmax, ymin), P2(xmax, ymax))
    SketchLine.Create(P2(xmax, ymax), P2(xmin, ymax))
    SketchLine.Create(P2(xmin, ymax), P2(xmin, ymin))

def sketch_doubleD_shifted(x0, y0, dx, R):
    """
    Double-D loop shifted by (x0,y0) in XY plane.
    (Same topology as your sketch_doubleD, but with offsets.)
    """
    # Right arc (center at +dx)
    SketchArc.CreateSweepArc(
        P2(x0 + dx, y0),
        P2(x0 + dx, y0 + R),
        P2(x0 + dx, y0 - R),
        True
    )

    # Left arc (center at -dx)
    SketchArc.CreateSweepArc(
        P2(x0 - dx, y0),
        P2(x0 - dx, y0 + R),
        P2(x0 - dx, y0 - R),
        False
    )

    # Top and bottom bridges
    SketchLine.Create(P2(x0 - dx, y0 + R), P2(x0 + dx, y0 + R))
    SketchLine.Create(P2(x0 - dx, y0 - R), P2(x0 + dx, y0 - R))

def extrude_largest_face_along_z(name, length_mm):
    solidify_sketch()

    # The "temp sketch body" is usually the last body created
    temp_body = GetRootPart().Bodies[-1]
    faces = list(temp_body.Faces)
    if not faces:
        raise Exception("No faces found to extrude for %s (sketch not closed?)" % name)

    # Pick largest face (for rectangle-minus-pocket, the largest should be the support region)
    face = max(faces, key=lambda f: f.Area)

    options = ExtrudeFaceOptions()
    options.ExtrudeType = ExtrudeType.ForceIndependent

    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirZ, MM(length_mm), options)
    created = list(res.CreatedBodies)
    if not created:
        raise Exception("Extrude created no bodies for %s" % name)

    body = created[0]
    body.SetName(name)
    return body

def move_body_to_component_and_refetch(body, target_component, final_name):
    """
    Move body into component and refetch by name (stable).
    Avoids relying on 'Bodies[-1]' which can be unstable.
    """
    # Ensure unique temp name before move
    body.SetName(final_name + "__TEMP__")

    sel = Selection.Create(body)
    ComponentHelper.MoveBodiesToComponent(sel, target_component)

    # Refetch by name inside the target component
    for b in list(target_component.Content.Bodies):
        try:
            if b.GetName() == final_name + "__TEMP__":
                b.SetName(final_name)
                return b
        except:
            if getattr(b, "Name", "") == final_name + "__TEMP__":
                b.SetName(final_name)
                return b

    raise Exception("Moved body but could not refetch '%s' in target component." % final_name)


def create_support_with_doubleD_pocket(name, z_center):
    # Outer overwrap envelope dims (double-D assumption)
    W_ow = W_outer + 2.0*t_overwrap
    H_ow = H_outer + 2.0*t_overwrap

    R_ow = 0.5 * H_ow
    dx_ow = 0.5 * (W_ow - H_ow)

    # Small clearance so cable isn't intersecting fixture
    clearance = 0.03 * min(W_ow, H_ow)
    R_pocket = R_ow + clearance
    dx_pocket = dx_ow

    # Support thickness along Z: derived + clamped
    t_support = 0.20 * L_extrude
    t_support = max(1.2*H_ow, min(t_support, 3.0*H_ow))

    # Block padding around envelope
    pad_x   = 0.35 * W_ow
    pad_y_up = 0.15 * H_ow
    pad_y_dn = 0.55 * H_ow

    # Pocket sits slightly below centerline so cable "rests"
    x0 = 0.0
    y0 = 0.10 * H_ow

    # Place by construction: start sketch at z0 so extrusion spans [z0, z0+t_support]
    z0 = z_center - 0.5*t_support
    set_sketch_plane_xy_at_z(z0)

    # Outer block rectangle
    xmin = -(0.5*W_ow + pad_x)
    xmax = +(0.5*W_ow + pad_x)
    ymin = -(0.5*H_ow + pad_y_dn)
    ymax = +(0.5*H_ow + pad_y_up)
    sketch_rectangle_xy(xmin, xmax, ymin, ymax)

    # Inner pocket loop
    sketch_doubleD_shifted(x0, y0, dx_pocket, R_pocket)

    temp_support = extrude_largest_face_along_z(name, t_support)
    return temp_support

# ============================================================
# 7B) Two Supports with Double-D Pocket (XY sketch -> Z extrude)
# ============================================================
def _body_name(b):
    try:
        return b.GetName()
    except:
        return getattr(b, "Name", "")

def find_bodies_by_prefix_anywhere(prefixes):
    """Return list of bodies in root + immediate components whose name starts with any prefix."""
    root = GetRootPart()
    out = []

    # root bodies
    for b in list(root.Bodies):
        n = _body_name(b)
        if any(n.startswith(p) for p in prefixes):
            out.append(b)

    # component bodies (1 level deep; enough for your Cable_Bodies / RigidParts_3Point_Bending)
    for comp in list(root.Components):
        try:
            for b in list(comp.Content.Bodies):
                n = _body_name(b)
                if any(n.startswith(p) for p in prefixes):
                    out.append(b)
        except:
            pass

    return out

def delete_bodies(bodies):
    if bodies:
        Delete.Execute(Selection.Create(bodies))
def _body_name(b):
    try:
        return b.GetName()
    except:
        return getattr(b, "Name", "")

def find_body_anywhere_by_name(name):
    root = GetRootPart()

    # root bodies
    for b in list(root.Bodies):
        if _body_name(b) == name:
            return b, None  # None means root

    # one level components (covers your Cable_Bodies / RigidParts_3Point_Bending)
    for comp in list(root.Components):
        try:
            for b in list(comp.Content.Bodies):
                if _body_name(b) == name:
                    return b, comp
        except:
            pass

    return None, None

def move_body_to_component_once(body, current_comp, target_comp, final_name):
    # already in the right place → do nothing
    if current_comp is target_comp:
        return body

    # if already named correctly, keep it; otherwise set a temp name for refetch
    tmp_name = final_name + "__TEMP__"
    body.SetName(tmp_name)

    ComponentHelper.MoveBodiesToComponent(Selection.Create(body), target_comp)

    # refetch moved body by temp name (small loop, only happens once)
    for b in list(target_comp.Content.Bodies):
        if _body_name(b) == tmp_name:
            b.SetName(final_name)
            return b

    # fallback: try final_name match
    for b in list(target_comp.Content.Bodies):
        if _body_name(b) == final_name:
            return b

    raise Exception("Move succeeded but could not refetch moved body: " + final_name)

def create_support_with_doubleD_pocket(name, z_center):
    # Outer overwrap envelope dims (double-D assumption)
    W_ow = W_outer + 2.0*t_overwrap
    H_ow = H_outer + 2.0*t_overwrap

    R_ow = 0.5 * H_ow
    dx_ow = 0.5 * (W_ow - H_ow)

    # Small clearance so cable isn't intersecting fixture
    clearance = 0.03 * min(W_ow, H_ow)
    R_pocket = R_ow + clearance
    dx_pocket = dx_ow

    # Support thickness along Z: derived + clamped
    t_support = 0.20 * L_extrude
    t_support = max(1.2*H_ow, min(t_support, 3.0*H_ow))

    # Block padding around envelope
    pad_x   = 0.35 * W_ow
    pad_y_up = 0.15 * H_ow
    pad_y_dn = 0.55 * H_ow

    # Pocket sits slightly below centerline so cable "rests"
    x0 = 0.0
    y0 = 0.0 #-0.10 * H_ow

    # Place by construction: start sketch at z0 so extrusion spans [z0, z0+t_support]
    z0 = z_center - 0.5*t_support
    set_sketch_plane_xy_at_z(z0)

    # Outer block rectangle
    xmin = -(0.5*W_ow + pad_x)
    xmax = +(0.5*W_ow + pad_x)
    ymin = -(0.5*H_ow + pad_y_dn)
    ymax = +(0.5*H_ow + pad_y_up)
    sketch_rectangle_xy(xmin, xmax, ymin, ymax)

    # Inner pocket loop
    sketch_doubleD_shifted(x0, y0, dx_pocket, R_pocket)

    temp_support = extrude_largest_face_along_z(name, t_support)
    return temp_support

def create_two_supports():
    z1 = 0.20 * L_extrude
    z2 = 0.80 * L_extrude

    rigid_comp = get_or_create_component("RigidParts_3Point_Bending")

    # --- If supports already exist, RETURN immediately (NO side-effects) ---
    b1, c1 = find_body_anywhere_by_name("Rig_Support_1")
    b2, c2 = find_body_anywhere_by_name("Rig_Support_2")

    if b1 is not None and b2 is not None:
        # Optional: ensure they live in the rigid component, but ONLY if needed
        b1 = move_body_to_component_once(b1, c1, rigid_comp, "Rig_Support_1")
        b2 = move_body_to_component_once(b2, c2, rigid_comp, "Rig_Support_2")

        # IMPORTANT: do NOT recreate named selections here during debugging
        return b1, b2

    # --- Create missing ones ONLY ---
    if b1 is None:
        s1_tmp = create_support_with_doubleD_pocket("Rig_Support_1", z1)
        b1 = move_body_to_component_once(s1_tmp, None, rigid_comp, "Rig_Support_1")
        # create_ns("Rig_Support_1", b1)  # add back later

    if b2 is None:
        s2_tmp = create_support_with_doubleD_pocket("Rig_Support_2", z2)
        b2 = move_body_to_component_once(s2_tmp, None, rigid_comp, "Rig_Support_2")
        # create_ns("Rig_Support_2", b2)  # add back later

    return b1, b2

support1, support2 = create_two_supports()

# ============================================================
# 7G) OPEN THE SUPPORTS: CUT WITH ZX PLANE (Y=0) AND KEEP BOTTOM
# ============================================================
# ============================================================
# 7G) OPEN THE SUPPORTS: SPLIT WITH ZX PLANE USING ByCutter(FACES)
#   Signature: ByCutter(bodySelection, toolFacesSelection, extendSurfaces, [info])
# ============================================================

def set_sketch_plane_zx_at_y(y0):
    # ZX plane at Y=y0. In-plane axes: Z (u), X (v)
    frame = Frame.Create(
        Point.Create(MM(0), MM(y0), MM(0)),
        Direction.DirZ,
        Direction.DirX
    )
    ViewHelper.SetSketchPlane(Plane.Create(frame))

def sketch_rectangle_uv(umin, umax, vmin, vmax):
    # Rectangle in current sketch plane coords (u,v) mapped to Point2D(x,y)
    SketchLine.Create(P2(umin, vmin), P2(umax, vmin))
    SketchLine.Create(P2(umax, vmin), P2(umax, vmax))
    SketchLine.Create(P2(umax, vmax), P2(umin, vmax))
    SketchLine.Create(P2(umin, vmax), P2(umin, vmin))

def make_thin_cutter_slab_ZX(y_cut=0.0, name="TMP_ZX_CUTTER"):
    """
    Make a very large, very thin slab centered at Y=y_cut.
    We'll use its large +/-Y faces as tool faces for SplitBody.ByCutter.
    """
    W_ow = W_outer + 2.0*t_overwrap
    H_ow = H_outer + 2.0*t_overwrap

    halfZ = 1.20 * L_extrude
    halfX = 1.50 * W_ow

    tY = max(0.02 * H_ow, 0.05)  # thin but nonzero

    # Build slab by extruding along +Y from plane at y_cut - tY/2
    set_sketch_plane_zx_at_y(y_cut - 0.5*tY)

    # In ZX sketch: u=Z, v=X
    sketch_rectangle_uv(-halfZ, +halfZ, -halfX, +halfX)

    solidify_sketch()
    temp_body = GetRootPart().Bodies[-1]
    face = max(list(temp_body.Faces), key=lambda f: f.Area)

    opts = ExtrudeFaceOptions()
    opts.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirY, MM(tY), opts)

    created = list(res.CreatedBodies)
    if not created:
        raise Exception("Failed to create cutter slab body (no bodies created).")

    slab = created[0]
    slab.SetName(name)
    return slab

def _faces_with_normal_near_dir(faces, dir_vec, tol=0.95):
    """
    Return faces whose normal aligns with dir_vec (dot > tol).
    Uses Face.Plane if available; falls back to picking by area if not.
    """
    out = []
    for f in faces:
        try:
            # Works for planar faces
            n = f.Plane.Normal
            # dot product
            d = n.X*dir_vec.X + n.Y*dir_vec.Y + n.Z*dir_vec.Z
            if d > tol:
                out.append(f)
        except:
            pass
    return out

def body_representative_y(body):
    """
    Returns a representative Y coordinate for the body by sampling a face centroid.
    Works across SpaceClaim builds (no GetBoundingBox needed).
    """
    for f in list(body.Faces):
        try:
            pt = f.Eval(0.5, 0.5).Point  # face parametric midpoint
            return pt.Y
        except:
            pass
    # Fallback: origin if something is very wrong
    return 0.0

def split_by_ZX_faces_keep_bottom(body, y_cut=0.0, final_name=None):
    """
    Split 'body' using SplitBody.ByCutter with tool *faces* for a ZX cut at Y=y_cut.
    Then delete ALL split pieces above, keeping exactly one (lowest representative Y).
    This version does NOT trust split_res.CreatedBodies; it instead:
      - renames the input body to a unique tag
      - after split, finds all bodies anywhere whose name starts with that tag
      - keeps the lowest-Y one, deletes the rest
      - deletes the cutter slab
    """

    # -------- helpers local to this function --------
    def body_representative_y(b):
        for f in list(b.Faces):
            try:
                pt = f.Eval(0.5, 0.5).Point
                return pt.Y
            except:
                pass
        return 0.0

    def _body_name(b):
        try:
            return b.GetName()
        except:
            return getattr(b, "Name", "")

    def find_bodies_by_prefix_anywhere(prefix):
        root = GetRootPart()
        out = []

        # root bodies
        for b in list(root.Bodies):
            if _body_name(b).startswith(prefix):
                out.append(b)

        # 1-level components
        for comp in list(root.Components):
            try:
                for b in list(comp.Content.Bodies):
                    if _body_name(b).startswith(prefix):
                        out.append(b)
            except:
                pass

        return out

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

    # -------- 1) Create thin slab cutter centered at y_cut --------
    slab = make_thin_cutter_slab_ZX(y_cut=y_cut, name="TMP_ZX_CUTTER")

    # -------- 2) Pick tool faces from slab (prefer +/-Y faces; fallback to two largest) --------
    slab_faces = list(slab.Faces)

    plusY = _faces_with_normal_near_dir(slab_faces, Direction.DirY)

    try:
        dir_minusY = Direction.DirY.Negate()
    except:
        dir_minusY = Direction.Create(0, -1, 0)

    minusY = _faces_with_normal_near_dir(slab_faces, dir_minusY)

    tool_faces = []
    if plusY:
        tool_faces.append(max(plusY, key=lambda f: f.Area))
    if minusY:
        tool_faces.append(max(minusY, key=lambda f: f.Area))

    if len(tool_faces) < 2:
        tool_faces = sorted(slab_faces, key=lambda f: f.Area, reverse=True)[:2]

    tool_faces_sel = Selection.Create(tool_faces)
    body_sel = Selection.Create(body)

    # -------- 3) Tag the body BEFORE split so we can find all split pieces reliably --------
    tag = (final_name if final_name else _body_name(body))
    if not tag:
        tag = "RigSupport"
    prefix = tag + "__SPLIT__"

    # Ensure uniqueness enough for this run
    body.SetName(prefix + "0")

    # -------- 4) Split using your signature --------
    extend_surfaces = True
    SplitBody.ByCutter(body_sel, tool_faces_sel, extend_surfaces)

    # -------- 5) Find all bodies produced from this split by name prefix --------
    pieces = find_bodies_by_prefix_anywhere(prefix)

    if len(pieces) < 2:
        # Split didn’t create multiple bodies (or naming behavior differs)
        # Cleanup cutter and return original body reference (which may be renamed)
        try:
            Delete.Execute(Selection.Create(slab))
        except:
            pass
        print("WARNING: Could not find split pieces by prefix; leaving body as-is.")
        if final_name:
            try:
                body.SetName(final_name)
            except:
                pass
        return body

    # -------- 6) Keep exactly ONE: the lowest-Y piece; delete all other pieces + slab --------
    keeper = None
    keeper_y = None
    for b in pieces:
        yb = body_representative_y(b)
        if keeper is None or yb < keeper_y:
            keeper = b
            keeper_y = yb

    to_delete = [b for b in pieces if b is not keeper]
    if to_delete:
        Delete.Execute(Selection.Create(to_delete))

    Delete.Execute(Selection.Create(slab))

    if final_name:
        keeper.SetName(final_name)
    else:
        keeper.SetName(tag)

    return keeper


# --- APPLY ---
support1 = split_by_ZX_faces_keep_bottom(support1, y_cut=0.0, final_name="Rig_Support_1")
support2 = split_by_ZX_faces_keep_bottom(support2, y_cut=0.0, final_name="Rig_Support_2")





################################################################################################
# ============================================================
# 7) RIGID 3-POINT BENDING: LOADING NOSE  +  TWO DOUBLE-D SUPPORTS
#   - Cable axis is Z (same as your extrusions)
#   - U/pocket profiles sketched in XY, extruded along Z
#   - Rigid bodies moved into component: RigidParts_3Point_Bending
#   - Supports located at 20% and 80% of L_extrude
# ============================================================

# -----------------------------
# Parameters (derived defaults)
# -----------------------------

# -----------------------------
# Helpers (minimal + stable)
# -----------------------------

def get_or_create_component(name):
    root = GetRootPart()
    # Find existing
    for c in list(root.Components):
        try:
            if c.GetName() == name:
                return c
        except:
            if getattr(c, "Name", "") == name:
                return c

    # Create via Part blueprint + Component instance (your working style)
    doc = Window.ActiveWindow.Document
    part_def = Part.Create(doc, name)
    comp = Component.Create(root, part_def)
    return comp

def move_body_to_component_and_refetch(body, target_component, final_name):
    """
    Moves 'body' into component and returns a stable reference by refetching via temp name.
    This avoids 'Bodies[-1]' instability and prevents ping-pong / regen weirdness.
    """
    tmp_name = final_name + "__TEMP__"
    body.SetName(tmp_name)

    ComponentHelper.MoveBodiesToComponent(Selection.Create(body), target_component)

    for b in list(target_component.Content.Bodies):
        try:
            if b.GetName() == tmp_name:
                b.SetName(final_name)
                return b
        except:
            if getattr(b, "Name", "") == tmp_name:
                b.SetName(final_name)
                return b

    raise Exception("Moved body but could not refetch '%s' in target component." % final_name)

def extrude_largest_face_along_y(name, length_mm):
    solidify_sketch()
    temp_body = GetRootPart().Bodies[-1]
    faces = list(temp_body.Faces)
    if not faces:
        raise Exception("No faces found to extrude for %s (sketch not closed?)" % name)
    face = max(faces, key=lambda f: f.Area)

    options = ExtrudeFaceOptions()
    options.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirY, MM(length_mm), options)

    created = list(res.CreatedBodies)
    if not created:
        raise Exception("Extrude created no bodies for %s" % name)

    body = created[0]
    body.SetName(name)
    return body



def sketch_rectangle_xy(xmin, xmax, ymin, ymax):
    SketchLine.Create(P2(xmin, ymin), P2(xmax, ymin))
    SketchLine.Create(P2(xmax, ymin), P2(xmax, ymax))
    SketchLine.Create(P2(xmax, ymax), P2(xmin, ymax))
    SketchLine.Create(P2(xmin, ymax), P2(xmin, ymin))

def sketch_doubleD_shifted(x0, y0, dx, R):
    # right arc
    SketchArc.CreateSweepArc(P2(x0 + dx, y0), P2(x0 + dx, y0 + R), P2(x0 + dx, y0 - R), True)
    # left arc
    SketchArc.CreateSweepArc(P2(x0 - dx, y0), P2(x0 - dx, y0 + R), P2(x0 - dx, y0 - R), False)
    # bridges
    SketchLine.Create(P2(x0 - dx, y0 + R), P2(x0 + dx, y0 + R))
    SketchLine.Create(P2(x0 - dx, y0 - R), P2(x0 + dx, y0 - R))

# ============================================================
# 7A) Loading Nose (YZ circle -> X extrude, centered by construction)
def move_body_to_component(body, target_component):
    """
    Moves a body to a component and RETURNS THE NEW BODY.
    Crucial: The original 'body' reference is invalid after the move.
    """
    sel = Selection.Create(body)
    ComponentHelper.MoveBodiesToComponent(sel, target_component)
    
    # The body is now the last one in the target component's list
    # We must return this NEW reference.
    return target_component.Content.Bodies[-1]

def extrude_last_profile_along_y(name, length_mm, pick_largest=True, dir_sign=+1):
    solidify_sketch()
    temp_body = GetRootPart().Bodies[-1]
    face = max(temp_body.Faces, key=lambda f: f.Area) if pick_largest else temp_body.Faces[0]

    options = ExtrudeFaceOptions()
    options.ExtrudeType = ExtrudeType.ForceIndependent

    diry = Direction.DirY if dir_sign > 0 else Direction.DirY.Negate()
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), diry, MM(length_mm), options)
    new_body = res.CreatedBodies[0]
    new_body.SetName(name)
    return new_body

def create_loading_nose_in_rigid_component_bad():
    r = Nose_Diam / 2.0

    # ... (Geometry Math is same) ...
    cable_top_x = 0.5 * (W_outer + 2.0 * t_overwrap)
    x_center = cable_top_x + Initial_Gap + r
    z_center = 0.5 * L_extrude

    # ... (Sketch & Extrude is same) ...
    set_sketch_plane_zx_at_y(-0.5 * Nose_Length)
    sketch_circle(x_center, z_center, r)
    
    # This 'nose' is temporary (lives in Root)
    temp_nose = extrude_last_profile_along_y("Rig_Loading_Nose_Bad", Nose_Length)

    # ... (Organize) ...
    rigid_comp = get_or_create_component("RigidParts_3Point_Bending_Bad")
    
    # FIX: Update 'nose' variable to the new body inside the component
    final_nose = move_body_to_component(temp_nose, rigid_comp)

    # ... (Named Selection) ...
    # Now we use 'final_nose', which is a valid reference
    create_ns("Rig_Loading_Nose_Bad", final_nose)

    return final_nose

# Execute
nose_body_bad = create_loading_nose_in_rigid_component_bad()

# ============================================================
# 7F) RIGID 3-POINT BENDING: TWO SUPPORTS WITH DOUBLE-D POCKET
#   - Pocket follows OUTER OVERWRAP shape (double-D envelope)
#   - Sketch in XY, extrude in Z (cable axis)
#   - Supports at 20% and 80% of L_extrude
#   - No boolean needed (pocket from sketch loops)
# ============================================================
def move_body_to_component_and_refetch(body, target_component, final_name):
    """
    Move body into component and refetch by name (stable).
    Avoids relying on 'Bodies[-1]' which can be unstable.
    """
    # Ensure unique temp name before move
    body.SetName(final_name + "__TEMP__")

    sel = Selection.Create(body)
    ComponentHelper.MoveBodiesToComponent(sel, target_component)

    # Refetch by name inside the target component
    for b in list(target_component.Content.Bodies):
        try:
            if b.GetName() == final_name + "__TEMP__":
                b.SetName(final_name)
                return b
        except:
            if getattr(b, "Name", "") == final_name + "__TEMP__":
                b.SetName(final_name)
                return b

    raise Exception("Moved body but could not refetch '%s' in target component." % final_name)

def set_sketch_plane_xy_at_z(z0):
    frame = Frame.Create(
        Point.Create(MM(0), MM(0), MM(z0)),
        Direction.DirX,
        Direction.DirY
    )
    plane = Plane.Create(frame)
    ViewHelper.SetSketchPlane(plane)

def sketch_rectangle_xy(xmin, xmax, ymin, ymax):
    SketchLine.Create(P2(xmin, ymin), P2(xmax, ymin))
    SketchLine.Create(P2(xmax, ymin), P2(xmax, ymax))
    SketchLine.Create(P2(xmax, ymax), P2(xmin, ymax))
    SketchLine.Create(P2(xmin, ymax), P2(xmin, ymin))


# ============================================================
# 7B) Two Supports with Double-D Pocket (XY sketch -> Z extrude)
# ============================================================


def create_two_supports_bad():
    z1 = 0.20 * L_extrude
    z2 = 0.80 * L_extrude

    rigid_comp = get_or_create_component("RigidParts_3Point_Bending_Bad")

    # --- If supports already exist, RETURN immediately (NO side-effects) ---
    b1, c1 = find_body_anywhere_by_name("Rig_Support_Bad_1")
    b2, c2 = find_body_anywhere_by_name("Rig_Support_Bad_2")

    if b1 is not None and b2 is not None:
        # Optional: ensure they live in the rigid component, but ONLY if needed
        b1 = move_body_to_component_once(b1, c1, rigid_comp, "Rig_Support_Bad_1")
        b2 = move_body_to_component_once(b2, c2, rigid_comp, "Rig_Support_Bad_2")

        # IMPORTANT: do NOT recreate named selections here during debugging
        return b1, b2

    # --- Create missing ones ONLY ---
    if b1 is None:
        s1_tmp = create_support_with_doubleD_pocket("Rig_Support_Bad_1", z1)
        b1 = move_body_to_component_once(s1_tmp, None, rigid_comp, "Rig_Support_Bad_1")
        # create_ns("Rig_Support_1", b1)  # add back later

    if b2 is None:
        s2_tmp = create_support_with_doubleD_pocket("Rig_Support_Bad_2", z2)
        b2 = move_body_to_component_once(s2_tmp, None, rigid_comp, "Rig_Support_Bad_2")
        # create_ns("Rig_Support_2", b2)  # add back later

    return b1, b2

support1_bad, support2_bad = create_two_supports_bad()

# ============================================================
# 7G) OPEN THE SUPPORTS: CUT WITH ZX PLANE (Y=0) AND KEEP BOTTOM
# ============================================================
# ============================================================
# 7G) OPEN THE SUPPORTS: SPLIT WITH ZX PLANE USING ByCutter(FACES)
#   Signature: ByCutter(bodySelection, toolFacesSelection, extendSurfaces, [info])
# ============================================================

def set_sketch_plane_zx_at_y(y0):
    # ZX plane at Y=y0. In-plane axes: Z (u), X (v)
    frame = Frame.Create(
        Point.Create(MM(0), MM(y0), MM(0)),
        Direction.DirZ,
        Direction.DirX
    )
    ViewHelper.SetSketchPlane(Plane.Create(frame))

def sketch_rectangle_uv(umin, umax, vmin, vmax):
    # Rectangle in current sketch plane coords (u,v) mapped to Point2D(x,y)
    SketchLine.Create(P2(umin, vmin), P2(umax, vmin))
    SketchLine.Create(P2(umax, vmin), P2(umax, vmax))
    SketchLine.Create(P2(umax, vmax), P2(umin, vmax))
    SketchLine.Create(P2(umin, vmax), P2(umin, vmin))

def make_thin_cutter_slab_YZ(x_cut=0.0, name="TMP_YZ_CUTTER"):
    """
    Make a very large, very thin slab centered at X=x_cut.
    We'll use its large +/-X faces as tool faces for SplitBody.ByCutter.
    """
    W_ow = W_outer + 2.0*t_overwrap
    H_ow = H_outer + 2.0*t_overwrap

    halfZ = 1.20 * L_extrude
    halfY = 1.50 * H_ow

    tX = max(0.02 * W_ow, 0.05)  # thin but nonzero

    # Build slab by extruding along +X from plane at x_cut - tX/2
    set_sketch_plane_yz_at_x(x_cut - 0.5*tX)

    # In YZ sketch: u=Z, v=Y
    sketch_rectangle_uv(-halfZ, +halfZ, -halfY, +halfY)

    solidify_sketch()
    temp_body = GetRootPart().Bodies[-1]
    face = max(list(temp_body.Faces), key=lambda f: f.Area)

    opts = ExtrudeFaceOptions()
    opts.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirX, MM(tX), opts)

    created = list(res.CreatedBodies)
    if not created:
        raise Exception("Failed to create cutter slab body (no bodies created).")

    slab = created[0]
    slab.SetName(name)
    return slab

def _faces_with_normal_near_dir(faces, dir_vec, tol=0.95):
    """
    Return faces whose normal aligns with dir_vec (dot > tol).
    Uses Face.Plane if available; falls back to picking by area if not.
    """
    out = []
    for f in faces:
        try:
            # Works for planar faces
            n = f.Plane.Normal
            # dot product
            d = n.X*dir_vec.X + n.Y*dir_vec.Y + n.Z*dir_vec.Z
            if d > tol:
                out.append(f)
        except:
            pass
    return out

def body_representative_x(body):
    """
    Returns a representative X coordinate for the body by sampling a face centroid.
    Works across SpaceClaim builds (no GetBoundingBox needed).
    """
    for f in list(body.Faces):
        try:
            pt = f.Eval(0.5, 0.5).Point  # face parametric midpoint
            return pt.X
        except:
            pass
    # Fallback: origin if something is very wrong
    return 0.0

def split_by_YZ_faces_keep_bottom(body, x_cut=0.0, final_name=None):
    """
    Split 'body' using SplitBody.ByCutter with tool *faces* for a YZ cut at x=x_cut.
    Then delete ALL split pieces above, keeping exactly one (lowest representative X).
    This version does NOT trust split_res.CreatedBodies; it instead:
      - renames the input body to a unique tag
      - after split, finds all bodies anywhere whose name starts with that tag
      - keeps the lowest-X one, deletes the rest
      - deletes the cutter slab
    """

    # -------- helpers local to this function --------
    def body_representative_x(b):
        for f in list(b.Faces):
            try:
                pt = f.Eval(0.5, 0.5).Point
                return pt.X
            except:
                pass
        return 0.0

    def _body_name(b):
        try:
            return b.GetName()
        except:
            return getattr(b, "Name", "")

    def find_bodies_by_prefix_anywhere(prefix):
        root = GetRootPart()
        out = []

        # root bodies
        for b in list(root.Bodies):
            if _body_name(b).startswith(prefix):
                out.append(b)

        # 1-level components
        for comp in list(root.Components):
            try:
                for b in list(comp.Content.Bodies):
                    if _body_name(b).startswith(prefix):
                        out.append(b)
            except:
                pass

        return out

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

    # -------- 1) Create thin slab cutter centered at x_cut --------
    slab = make_thin_cutter_slab_YZ(x_cut=x_cut, name="TMP_YZ_CUTTER")

    # -------- 2) Pick tool faces from slab (prefer +/-X faces; fallback to two largest) --------
    slab_faces = list(slab.Faces)

    plusX = _faces_with_normal_near_dir(slab_faces, Direction.DirX)

    try:
        dir_minusX = Direction.DirX.Negate()
    except:
        dir_minusX = Direction.Create(-1, 0, 0)

    minusX = _faces_with_normal_near_dir(slab_faces, dir_minusX)

    tool_faces = []
    if plusX:
        tool_faces.append(max(plusX, key=lambda f: f.Area))
    if minusX:
        tool_faces.append(max(minusX, key=lambda f: f.Area))

    if len(tool_faces) < 2:
        tool_faces = sorted(slab_faces, key=lambda f: f.Area, reverse=True)[:2]

    tool_faces_sel = Selection.Create(tool_faces)
    body_sel = Selection.Create(body)

    # -------- 3) Tag the body BEFORE split so we can find all split pieces reliably --------
    tag = (final_name if final_name else _body_name(body))
    if not tag:
        tag = "RigSupport"
    prefix = tag + "__SPLIT__"

    # Ensure uniqueness enough for this run
    body.SetName(prefix + "0")

    # -------- 4) Split using your signature --------
    extend_surfaces = True
    SplitBody.ByCutter(body_sel, tool_faces_sel, extend_surfaces)

    # -------- 5) Find all bodies produced from this split by name prefix --------
    pieces = find_bodies_by_prefix_anywhere(prefix)

    if len(pieces) < 2:
        # Split didn’t create multiple bodies (or naming behavior differs)
        # Cleanup cutter and return original body reference (which may be renamed)
        try:
            Delete.Execute(Selection.Create(slab))
        except:
            pass
        print("WARNING: Could not find split pieces by prefix; leaving body as-is.")
        if final_name:
            try:
                body.SetName(final_name)
            except:
                pass
        return body

    # -------- 6) Keep exactly ONE: the lowest-X piece; delete all other pieces + slab --------
    keeper = None
    keeper_x = None
    for b in pieces:
        xb = body_representative_x(b)
        if keeper is None or xb < keeper_x:
            keeper = b
            keeper_x = xb

    to_delete = [b for b in pieces if b is not keeper]
    if to_delete:
        Delete.Execute(Selection.Create(to_delete))

    Delete.Execute(Selection.Create(slab))

    if final_name:
        keeper.SetName(final_name)
    else:
        keeper.SetName(tag)

    return keeper


# --- APPLY ---
support1 = split_by_YZ_faces_keep_bottom(support1_bad, x_cut=0.0, final_name="Rig_Support_Bad_1")
support2 = split_by_YZ_faces_keep_bottom(support2_bad, x_cut=0.0, final_name="Rig_Support_Bad_2")



