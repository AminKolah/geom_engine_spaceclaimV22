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
# 7) RIGID 3-POINT BENDING: LOADING NOSE (Fixed)
# ============================================================
Nose_Diam = get_param("Nose_Diam", 2 * H_outer)
Nose_Length = get_param("Nose_Length", 4 * H_outer)
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
# 8) RIGID 3-POINT BENDING: SUPPORTS (Left & Right)
# ============================================================

# Parameters for Supports
Support_Diam = get_param("Support_Diam", H_outer)
Support_Span = get_param("Support_Span", 2 * H_outer) # Distance between supports
# Re-using Initial_Gap from Step 7, or you can define a new one
Support_Gap  = get_param("Support_Gap", 0.5) 

def create_supports_in_rigid_component():
    r = Support_Diam / 2.0
    
    # ... (Geometry calculations remain the same) ...
    cable_bot_y = -0.5 * (H_outer + 2.0 * t_overwrap)
    y_center = cable_bot_y - Support_Gap - r
    z_mid = 0.5 * L_extrude
    z_left  = z_mid - (0.5 * Support_Span)
    z_right = z_mid + (0.5 * Support_Span)
    
    rigid_comp = get_or_create_component("RigidParts_3Point_Bending")
    
    created_supports = []
    
    definitions = [
        ("Rig_Support_Left",  z_left),
        ("Rig_Support_Right", z_right)
    ]
    
    for name, z_pos in definitions:
        # FIX: Set the sketch plane INSIDE the loop.
        # This forces SpaceClaim back into Sketch Mode for every new support.
        set_sketch_plane_yz_at_x(-0.5 * Nose_Length) 
        
        # 1. Sketch
        sketch_circle(y_center, z_pos, r)
        
        # 2. Extrude (This switches us to 3D mode)
        temp_body = extrude_last_profile_along_x(name, Nose_Length)
        
        # 3. Move & Update Reference
        final_body = move_body_to_component(temp_body, rigid_comp)
        
        # 4. Name Selection
        create_ns(name, final_body)
        
        created_supports.append(final_body)

    return created_supports
# Execute
supports = create_supports_in_rigid_component()