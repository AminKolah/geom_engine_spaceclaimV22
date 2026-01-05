# Python Script, API Version = V22

from SpaceClaim.Api.V22 import *
from SpaceClaim.Api.V22.Geometry import *

# -----------------------------
# EXPECTED SCRIPT PARAMETERS (define in Groups -> Script Parameters):
#   R_inner  (Length, mm)
#   t_shell  (Length, mm)
# OPTIONAL:
#   L        (Length, mm)  [otherwise uses default below]
# -----------------------------

# Default length if you didn't create L as a script parameter
try:
    L_val = float(Parameters.L)   # if Script Parameter "L" exists
except:
    L_val = 7.0

# Get radii from Script Parameters
try:
    R_inner_val = float(Parameters.R_inner)
    t_shell_val = float(Parameters.t_shell)
except NameError:
    raise Exception(
        "Script parameters not found. Create Script Parameters named "
        "'R_inner' and 't_shell' under Groups -> Script Parameters (Length, mm)."
    )

R_outer_val = R_inner_val + t_shell_val

if R_inner_val <= 0:
    raise Exception("R_inner must be > 0.")
if t_shell_val <= 0:
    raise Exception("t_shell must be > 0.")

# -----------------------------
# Helpers
# -----------------------------
def set_sketch_plane_xy():
    ViewHelper.SetSketchPlane(Plane.PlaneXY)

def solidify_sketch():
    ViewHelper.SetViewMode(InteractionMode.Solid)

def sketch_circle_origin(r_mm):
    origin = Point2D.Create(MM(0), MM(0))
    SketchCircle.Create(origin, MM(r_mm))
    
def extrude_face(target_face, length, independent=True):
    # Create the required selection object
    selection = FaceSelection.Create(target_face)
    
    # Initialize extrusion options
    options = ExtrudeFaceOptions()
    
    # Logic to switch between Merging and No Merge
    if independent:
        options.ExtrudeType = ExtrudeType.ForceIndependent
    else:
        options.ExtrudeType = ExtrudeType.Add # Standard behavior
    
    # Execute the command
    return ExtrudeFaces.Execute(selection, length, options)

#def extrude_face(face_obj, length_mm):
 #   options = ExtrudeFaceOptions()
#    options.ExtrudeType = ExtrudeType.ForceIndependent
 #   ExtrudeFaces.Execute(face_obj, MM(length_mm), options)

# ==========================================================
# (1) INNER CYLINDER
# ==========================================================
set_sketch_plane_xy()
sketch_circle_origin(R_inner_val)
result = solidify_sketch()
FACE_INNER = GetRootPart().Bodies[0].Faces[0]
extrude_face(FACE_INNER, L_val)

# ==========================================================
# (2) OUTER TUBE (annulus)
# ==========================================================
set_sketch_plane_xy()
sketch_circle_origin(R_outer_val)
sketch_circle_origin(R_inner_val)
result = solidify_sketch()
FACE_ANNULUS = GetRootPart().Bodies[1].Faces[0]
extrude_face(FACE_ANNULUS, L_val)

print(
    "Done: inner radius =",
    R_inner_val,
    " thickness =",
    t_shell_val,
    " outer radius =",
    R_outer_val,
    " length =",
    L_val
)



# Rename 'Solid' to 'Inner_Cylinder'
body = GetRootPart().Bodies[0]
# 2. Create a Selection object from that body
selection = Selection.Create(body)
# 3. Rename using the Execute command
RenameObject.Execute(selection, "Inner_Cylinder")

# Rename 'Solid' to 'Inner_Cylinder'
body = GetRootPart().Bodies[1]
# 2. Create a Selection object from that body
selection = Selection.Create(body)
# 3. Rename using the Execute command
RenameObject.Execute(selection, "Outer_shell")
