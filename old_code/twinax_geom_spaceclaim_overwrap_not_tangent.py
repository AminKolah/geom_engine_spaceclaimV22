# Python Script for SpaceClaim V22
# Robust Cable Geometry Generation with Drains and Overwrap
from SpaceClaim.Api.V22 import *
from SpaceClaim.Api.V22.Geometry import *
from SpaceClaim.Api.V22.Modeler import *

# -----------------------------
# 1. PARAMETERS / INPUTS
# -----------------------------
def get_param(name, default_val):
    try: return float(Parameters[name])
    except: return default_val

# Conductor and Core
D_cond    = get_param('D_cond', 1.0)
C2C       = get_param('C2C', 3.0)
D_core    = get_param('D_core', 3.0)

# Shield (Double-D)
t_shield  = get_param('t_shield', 0.1)
W_outer   = get_param('W_outer', 8.0)
H_outer   = get_param('H_outer', 5.0)

# Drains and Overwrap
D_drain   = get_param('D_drain', 1.0)
t_overwrap = get_param('t_overwrap', 0.3)
L_extrude = get_param('L_extrude', 50.0)

# Derived values
r_cond = D_cond / 2.0
r_core = D_core / 2.0
r_drain = D_drain / 2.0
dx_cond = C2C / 2.0

# Shield inner dimensions
W_in = W_outer - 2.0 * t_shield
H_in = H_outer - 2.0 * t_shield
R_in = H_in / 2.0
dx_in = (W_in - H_in) / 2.0
R_shield = H_outer / 2.0
dx_shield = (W_outer - H_outer) / 2.0

ClearAll()

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
set_sketch_plane_xy()

# (1) Conductors
sketch_circle(-dx_cond, 0, r_cond)
extrude_and_name("conductor[1]", L_extrude)

set_sketch_plane_xy()
sketch_circle(dx_cond, 0, r_cond)
extrude_and_name("conductor[2]", L_extrude)

# (2) Cores (Insulation)
set_sketch_plane_xy()
sketch_circle(-dx_cond, 0, r_core)
sketch_circle(-dx_cond, 0, r_cond)
extrude_and_name("single_core[1]", L_extrude, True)

set_sketch_plane_xy()
sketch_circle(dx_cond, 0, r_core)
sketch_circle(dx_cond, 0, r_cond)
extrude_and_name("single_core[2]", L_extrude, True)

# (3) Filler (Inner Space)
set_sketch_plane_xy()
sketch_doubleD(dx_in, R_in)
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
set_sketch_plane_xy()
sketch_circle(-x_drain, 0, r_drain)
extrude_and_name("drain[1]", L_extrude)

set_sketch_plane_xy()
sketch_circle(x_drain, 0, r_drain)
extrude_and_name("drain[2]", L_extrude)

# (6) Overwrap
# Scale overwrap width to accommodate the shield + drains
W_over_in = (x_drain + r_drain) * 2
dx_over_in = (W_over_in - H_outer) / 2.0
dx_over_out = dx_over_in 
R_over_out = (H_outer / 2.0) + t_overwrap

set_sketch_plane_xy()
sketch_doubleD(dx_over_out, R_over_out)
sketch_doubleD(dx_over_in, H_outer / 2.0)
extrude_and_name("Overwrap", L_extrude, True)

# -----------------------------
# 4. FINALIZE FOR ANSYS
# -----------------------------
for body in GetRootPart().Bodies:
    body.SharedTopology = SharedTopologyType.Share

print("Full Cable Assembly Created Successfully.")