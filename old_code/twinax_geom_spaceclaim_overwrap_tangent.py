# Python Script, API Version = V22
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
"""
W_over_in = (x_drain + r_drain) * 2
dx_over_in = (W_over_in - H_outer) / 2.0
dx_over_out = dx_over_in 
R_over_out = (H_outer / 2.0) + t_overwrap

set_sketch_plane_xy()
sketch_doubleD(dx_over_out, R_over_out)
sketch_doubleD(dx_over_in, H_outer / 2.0)
extrude_and_name("Overwrap", L_extrude, True)
"""
# ---------------------------------------------------------
# (9) REALISTIC TANGENT OVERWRAP
# ---------------------------------------------------------
set_sketch_plane_xy()

# Geometry Data
# R_shield (Large radius), r_drain (Small radius)
# distance 'd' between center of drain and center of shield arc
d = abs(x_drain - dx_shield) 
R = R_shield
r = r_drain

# Angle of the tangent line relative to the horizontal
# Calculation: sin(theta) = (R - r) / d
theta = math.asin((R - r) / d)

# Tangent Points on Shield Arc (relative to shield arc center)
tx_S = R * math.sin(theta)
ty_S = R * math.cos(theta)

# Tangent Points on Drain Arc (relative to drain center)
tx_D = r * math.sin(theta)
ty_D = r * math.cos(theta)

def sketch_wrap_boundary(offset):
    # Adjust radii for thickness/offset
    Ro = R + offset
    ro = r + offset
    
    # Recalculate theta for the offset (though usually same)
    # x_drain_abs is the drain center
    
    # 1. Arcs around the drains (outer side)
    # From top tangent point to bottom tangent point
    # Right Drain
    SketchArc.CreateSweepArc(P2(x_drain, 0), P2(x_drain + tx_D, ty_D), P2(x_drain + tx_D, -ty_D), True)
    # Left Drain
    SketchArc.CreateSweepArc(P2(-x_drain, 0), P2(-x_drain - tx_D, ty_D), P2(-x_drain - tx_D, -ty_D), False)
    
    # 2. Arcs on the Shield sides (inner side of wrap)
    # Right Shield Arc
    SketchArc.CreateSweepArc(P2(dx_shield, 0), P2(dx_shield + tx_S, ty_S), P2(dx_shield + tx_S, -ty_S), True)
    # Left Shield Arc
    SketchArc.CreateSweepArc(P2(-dx_shield, 0), P2(-dx_shield - tx_S, ty_S), P2(-dx_shield - tx_S, -ty_S), False)
    
    # 3. Connecting Tangent Lines (The 4 corners)
    # Top Right
    SketchLine.Create(P2(x_drain + tx_D, ty_D), P2(dx_shield + tx_S, ty_S))
    # Bottom Right
    SketchLine.Create(P2(x_drain + tx_D, -ty_D), P2(dx_shield + tx_S, -ty_S))
    # Top Left
    SketchLine.Create(P2(-x_drain - tx_D, ty_D), P2(-dx_shield - tx_S, ty_S))
    # Bottom Left
    SketchLine.Create(P2(-x_drain - tx_D, -ty_D), P2(-dx_shield - tx_S, -ty_S))

    # 4. Top and Bottom Flats (Connecting the two shield arcs)
    SketchLine.Create(P2(-dx_shield - tx_S, ty_S), P2(dx_shield + tx_S, ty_S))
    SketchLine.Create(P2(-dx_shield - tx_S, -ty_S), P2(dx_shield + tx_S, -ty_S))

# Draw inner boundary (0 offset) and outer boundary (t_overwrap offset)
sketch_wrap_boundary(0)
sketch_wrap_boundary(t_overwrap)

# Extrude
extrude_and_name("Overwrap", L_extrude, True)
# -----------------------------
# 4. FINALIZE FOR ANSYS
# -----------------------------
#for body in GetRootPart().Bodies:
#    body.SharedTopology = SharedTopologyType.Share

print("Full Cable Assembly Created Successfully.")