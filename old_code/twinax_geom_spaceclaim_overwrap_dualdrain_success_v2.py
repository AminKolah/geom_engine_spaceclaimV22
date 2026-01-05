# Python Script, API Version = V22
# Python Script for SpaceClaim V22
# Robust Cable Geometry Generation with or without Drains, shell or fill options
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
# filler_opt: 0 = fill, 1 = shell
# drain_opt:  0 = no drains, 1 = with drains
filler_opt  = int(get_param('filler_option', 1))
drain_opt   = int(get_param('has_drains', 0))
is_a_doublet = int(get_param('is_a_doublet', 1)) # 0 = Separate, 1 = Combined

# filler_mode logic
filler_mode = "shell" if filler_opt == 1 else "fill"

# Conductor and Core
D_cond    = get_param('diam_cond', 1.0)
C2C       = get_param('c2c', 3.0)
D_core    = get_param('diam_core', 3.0)

# Shield (Double-D)
t_shield  = get_param('t_shield', 0.1)
W_outer   = get_param('W_outer', 8.0)
H_outer   = get_param('H_outer', 5.0)

# Drains and Overwrap
D_drain   = get_param('diam_drain', 1.0)
t_overwrap = get_param('t_overwrap', 0.1)
L_extrude = get_param('length_extrude', 50.0)

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
filler_mode = "shell" # Options: "fill" or "shell"


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
    # Right semicircle
    c1 = SketchArc.CreateSweepArc(Point2D.Create(MM(dx), 0), Point2D.Create(MM(dx), MM(R)), Point2D.Create(MM(dx), MM(-R)), True).CreatedCurves
    # Left semicircle
    c2 = SketchArc.CreateSweepArc(Point2D.Create(MM(-dx), 0), Point2D.Create(MM(-dx), MM(R)), Point2D.Create(MM(-dx), MM(-R)), False).CreatedCurves
    # Connectors
    l1 = SketchLine.Create(Point2D.Create(MM(-dx), MM(R)), Point2D.Create(MM(dx), MM(R))).CreatedCurves
    l2 = SketchLine.Create(Point2D.Create(MM(-dx), MM(-R)), Point2D.Create(MM(dx), MM(-R))).CreatedCurves
    return list(c1) + list(c2) + list(l1) + list(l2)

def extrude_and_name(name, length_mm, pick_largest=False):
    solidify_sketch()
    temp_body = GetRootPart().Bodies[-1]
    
    # Select face (largest for annulus/filler, first for simple circles)
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

# (1) Conductors - Always Created
set_sketch_plane_xy()
sketch_circle(-dx_cond, 0, r_cond)
extrude_and_name("conductor[1]", L_extrude)

set_sketch_plane_xy()
sketch_circle(dx_cond, 0, r_cond)
extrude_and_name("conductor[2]", L_extrude)

# (2) Cores (Insulation) - Only if NOT a doublet
if not (is_a_doublet):
    set_sketch_plane_xy()
    sketch_circle(-dx_cond, 0, r_core)
    sketch_circle(-dx_cond, 0, r_cond)
    extrude_and_name("single_core[1]", L_extrude, True)

    set_sketch_plane_xy()
    sketch_circle(dx_cond, 0, r_core)
    sketch_circle(dx_cond, 0, r_cond)
    extrude_and_name("single_core[2]", L_extrude, True)

# (3) Filler / Doublet (Inner Space)
set_sketch_plane_xy()

if is_a_doublet and filler_opt == 1:
    # DOUBLET MODE: Combine Core and Shell into one body
    sketch_doubleD(dx_in, R_in)
    sketch_circle(-dx_cond, 0, r_cond) # Subtract only conductor
    sketch_circle(dx_cond, 0, r_cond)
    extrude_and_name("Second_Extrusion", L_extrude, True)

elif filler_mode.lower() == "fill":
    # Standard: Shield Inner minus Core Holes
    sketch_doubleD(dx_in, R_in)
    sketch_circle(-dx_cond, 0, r_core)
    sketch_circle(dx_cond, 0, r_core)
    extrude_and_name("Second_Extrusion", L_extrude, True)

elif filler_mode.lower() == "shell":
    # Shell: A band between the shield and an inner offset
    t_filler_shell = (dx_in + R_in) - (dx_cond + r_core)
    sketch_doubleD(dx_in, R_in)
    R_fill_inner = R_in - t_filler_shell
    if R_fill_inner > 0:
        sketch_doubleD(dx_in, R_fill_inner)
    
    sketch_circle(-dx_cond, 0, r_core)
    sketch_circle(dx_cond, 0, r_core)
    extrude_and_name("Second_Extrusion", L_extrude, True)

# (4) Shield
set_sketch_plane_xy()
sketch_doubleD(dx_shield, R_shield)
sketch_doubleD(dx_in, R_in)
extrude_and_name("Shield", L_extrude, True)

# (5) Drains (Tangent to Shield)
x_drain = dx_shield + R_shield + r_drain
if drain_opt:
    set_sketch_plane_xy()
    sketch_circle(-x_drain, 0, r_drain)
    extrude_and_name("drain[1]", L_extrude)
    set_sketch_plane_xy()
    sketch_circle(x_drain, 0, r_drain)
    extrude_and_name("drain[2]", L_extrude)

# (6) Overwrap

def P2(x, y):
    return Point2D.Create(MM(x), MM(y))

def sketch_wrap_loop_with_drains(offset):
    # Core Geometry for Overwrap
    d = abs(x_drain - dx_shield) 
    R = R_shield
    r = r_drain
    # The tangent angle is the same for any concentric offset
    theta = math.asin((R - r) / d)
    # Effective radii with offset
    Ro = R + offset
    ro = r + offset
    
    # Calculate Tangent Points for this specific offset
    # Shield side
    tx_S = Ro * math.sin(theta)
    ty_S = Ro * math.cos(theta)
    # Drain side
    tx_D = ro * math.sin(theta)
    ty_D = ro * math.cos(theta)

    # 1. Arcs around the Drains (Outer side)
    # Right: From top tangent to bottom tangent
    SketchArc.CreateSweepArc(P2(x_drain, 0), P2(x_drain + tx_D, ty_D), P2(x_drain + tx_D, -ty_D), True)
    # Left: From top tangent to bottom tangent
    SketchArc.CreateSweepArc(P2(-x_drain, 0), P2(-x_drain - tx_D, ty_D), P2(-x_drain - tx_D, -ty_D), False)

    # 2. Tangent Lines (The "stretching" parts)
    SketchLine.Create(P2(x_drain + tx_D, ty_D), P2(dx_shield + tx_S, ty_S))      # Top Right
    SketchLine.Create(P2(x_drain + tx_D, -ty_D), P2(dx_shield + tx_S, -ty_S))    # Bottom Right
    SketchLine.Create(P2(-x_drain - tx_D, ty_D), P2(-dx_shield - tx_S, ty_S))    # Top Left
    SketchLine.Create(P2(-x_drain - tx_D, -ty_D), P2(-dx_shield - tx_S, -ty_S))  # Bottom Left

    # 1. Top-Right Shoulder (Connects Tangent to Top Flat)
    # From tangent point (~2 o'clock) to top center (~12 o'clock)
    SketchArc.CreateSweepArc(P2(dx_shield, 0), P2(dx_shield + tx_S, ty_S), P2(dx_shield, Ro), False)

    # 2. Bottom-Right Shoulder (Connects Tangent to Bottom Flat)
    # From tangent point (~4 o'clock) to bottom center (~6 o'clock)
    SketchArc.CreateSweepArc(P2(dx_shield, 0), P2(dx_shield + tx_S, -ty_S), P2(dx_shield, -Ro), True)

    # 3. Top-Left Shoulder
    # From top center (~12 o'clock) to tangent point (~10 o'clock)
    SketchArc.CreateSweepArc(P2(-dx_shield, 0), P2(-dx_shield, Ro), P2(-dx_shield - tx_S, ty_S), False)

    # 4. Bottom-Left Shoulder
    # From bottom center (~6 o'clock) to tangent point (~8 o'clock)
    SketchArc.CreateSweepArc(P2(-dx_shield, 0), P2(-dx_shield, -Ro), P2(-dx_shield - tx_S, -ty_S), True)
    # 4. Top and Bottom Flats (Sitting exactly at Y = +/- (H_outer/2 + offset))
    
    SketchLine.Create(P2(-dx_shield, Ro), P2(dx_shield, Ro))
    SketchLine.Create(P2(-dx_shield, -Ro), P2(dx_shield, -Ro))

if drain_opt:
    set_sketch_plane_xy()
    sketch_wrap_loop_with_drains(0)             # Hugging the cable
    sketch_wrap_loop_with_drains(t_overwrap)    # Outer skin
    extrude_and_name("Overwrap", L_extrude, True)
else:
    set_sketch_plane_xy()
    # 1. Sketch the OUTER profile
    sketch_doubleD(dx_shield, R_shield + t_overwrap)
    # 2. Sketch the INNER profile (to create the hollow annulus)
    sketch_doubleD(dx_shield, R_shield)
    extrude_and_name("Overwrap", L_extrude, True)


# -----------------------------
# 4. FINALIZE & CLEANUP
# -----------------------------


#for body in GetRootPart().Bodies:
 #   body.SharedTopology = SharedTopologyType.Share

print("Full Cable Assembly Created Successfully. Doublet mode: " + str(bool(is_a_doublet)))