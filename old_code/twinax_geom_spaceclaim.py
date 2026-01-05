# Python Script, API Version = V22
# Python Script, API Version = V22
from SpaceClaim.Api.V22 import *
from SpaceClaim.Api.V22.Geometry import *

# -----------------------------
# INPUTS (mm)
# -----------------------------
D_cond    = 1.0
C2C       = 3.0

# core OD (dielectric around conductor)
core_mode = "tangent"    # "tangent" or "explicit"
D_core    = C2C          # used in explicit mode
gap_core  = 0.0          # used in tangent mode (0 => tangent)

# Outer shield envelope (outer boundary)
t_shield  = 0.10
W_outer   = 8.0
H_outer   = 5.0

# Filler behavior
filler_mode = "fill"     # "fill" or "shell"
t_shell     = 0.5        # only used for filler_mode="shell" (shell thickness inside shield inner boundary)

# Extrusion length
L_extrude = 200.0

# -----------------------------
# Helpers
# -----------------------------
def set_sketch_plane_xy():
    ViewHelper.SetSketchPlane(Plane.PlaneXY)

def solidify_sketch():
    ViewHelper.SetViewMode(InteractionMode.Solid)

def P2(x, y):
    return Point2D.Create(MM(x), MM(y))

def sketch_circle(cx, cy, r):
    SketchCircle.Create(P2(cx, cy), MM(r))

def sketch_doubleD(dx, R):
    # Right semicircle
    origin = P2(+dx, 0.0)
    start  = P2(+dx, +R)
    end    = P2(+dx, -R)
    SketchArc.CreateSweepArc(origin, start, end, True)

    # Left semicircle
    origin = P2(-dx, 0.0)
    start  = P2(-dx, +R)
    end    = P2(-dx, -R)
    SketchArc.CreateSweepArc(origin, start, end, False)

    # Connectors
    SketchLine.Create(P2(-dx, +R), P2(+dx, +R))
    SketchLine.Create(P2(-dx, -R), P2(+dx, -R))

def extrude_face_from_last_solidified_region(face_pick_index, length_mm, independent=True):
    """
    After Solidify, SpaceClaim typically creates a temporary planar body with faces for the regions.
    We grab the *last* body and select one of its faces by index.
    If you extrude the wrong region, change face_pick_index.
    """
    b = GetRootPart().Bodies[-1]
    faces = list(b.Faces)
    if len(faces) == 0:
        raise Exception("No faces found after Solidify. Region creation failed.")

    if face_pick_index < 0 or face_pick_index >= len(faces):
        raise Exception("face_pick_index out of range. Have %d faces; requested %d" % (len(faces), face_pick_index))

    target_face = faces[face_pick_index]
    selection = FaceSelection.Create(target_face)

    options = ExtrudeFaceOptions()
    options.ExtrudeType = ExtrudeType.ForceIndependent if independent else ExtrudeType.Add

    return ExtrudeFaces.Execute(selection, MM(length_mm), options)

def sketch_and_extrude(face_pick_index, length_mm):
    solidify_sketch()
    extrude_face_from_last_solidified_region(face_pick_index, length_mm, independent=True)

# -----------------------------
# Derived geometry
# -----------------------------
r_cond = D_cond / 2.0

# core radius
if core_mode.lower() == "tangent":
    D_core_eff = C2C + gap_core
else:
    D_core_eff = D_core

r_core = D_core_eff / 2.0
dx_cond = C2C / 2.0

# Outer shield double-D
if W_outer < H_outer:
    raise Exception("Require W_outer >= H_outer for double-D shape.")
R_shield  = H_outer / 2.0
dx_shield = (W_outer - H_outer) / 2.0

# Inner shield boundary (also outer boundary for filler)
W_in = W_outer - 2.0*t_shield
H_in = H_outer - 2.0*t_shield
if W_in <= 0 or H_in <= 0:
    raise Exception("t_shield too large.")
R_in  = H_in / 2.0
dx_in = (W_in - H_in) / 2.0

# For filler shell mode: inner boundary for shell thickness
W_shell_in = W_in - 2.0*t_shell
H_shell_in = H_in - 2.0*t_shell
if filler_mode.lower() == "shell":
    if W_shell_in <= 0 or H_shell_in <= 0:
        raise Exception("t_shell too large for inner shield boundary.")
R_shell_in  = H_shell_in / 2.0
dx_shell_in = (W_shell_in - H_shell_in) / 2.0

# -----------------------------
# BUILD GEOMETRY (sequential sketches)
# -----------------------------
set_sketch_plane_xy()

# (1) Left conductor (disk -> cylinder)
sketch_circle(-dx_cond, 0.0, r_cond)
# Usually the disk is Face1 (or Face0). In your clean examples it was Face1.
# Here we use index 0 to be safer when only one region exists.
sketch_and_extrude(face_pick_index=0, length_mm=L_extrude)

# (2) Left core (annulus -> tube)
set_sketch_plane_xy()
sketch_circle(-dx_cond, 0.0, r_core)
sketch_circle(-dx_cond, 0.0, r_cond)
# Annulus is usually a different face than the inner disk.
# In a clean doc, often Face0=annulus, Face1=disk (or vice versa).
# Start with 0; if it extrudes the wrong thing, change to 1.
sketch_and_extrude(face_pick_index=0, length_mm=L_extrude)

# (3) Right conductor
set_sketch_plane_xy()
sketch_circle(+dx_cond, 0.0, r_cond)
sketch_and_extrude(face_pick_index=0, length_mm=L_extrude)

# (4) Right core (annulus)
set_sketch_plane_xy()
sketch_circle(+dx_cond, 0.0, r_core)
sketch_circle(+dx_cond, 0.0, r_cond)
sketch_and_extrude(face_pick_index=0, length_mm=L_extrude)

# (5) Filler
# Trick to avoid booleans:
#   Sketch the filler OUTER boundary, sketch the core OUTER circles as holes,
#   then extrude the *single face* corresponding to (outer minus holes).
set_sketch_plane_xy()

mode = filler_mode.lower()
if mode == "fill":
    sketch_doubleD(dx_in, R_in)                 # filler outer boundary (inside shield)
    sketch_circle(-dx_cond, 0.0, r_core)        # holes
    sketch_circle(+dx_cond, 0.0, r_core)
    # The "outer-with-holes" region is typically the largest face. Often it's index 0, but not guaranteed.
    # Start with 0; if you extrude the wrong region, try 1 or 2.
    sketch_and_extrude(face_pick_index=0, length_mm=L_extrude)

elif mode == "shell":
    # shell = (inner shield boundary ring), not filling the center region
    # We sketch outer = inner shield boundary; inner = shrunken boundary; plus core holes
    sketch_doubleD(dx_in, R_in)                 # outer
    sketch_doubleD(dx_shell_in, R_shell_in)     # inner
    sketch_circle(-dx_cond, 0.0, r_core)        # holes
    sketch_circle(+dx_cond, 0.0, r_core)
    sketch_and_extrude(face_pick_index=0, length_mm=L_extrude)

else:
    raise Exception("Unknown filler_mode: %s (use 'fill' or 'shell')" % filler_mode)

# (6) Shield (double-D ring)
# Same no-boolean trick: outer boundary + inner boundary => ring face
set_sketch_plane_xy()
sketch_doubleD(dx_shield, R_shield)  # outer shield boundary
sketch_doubleD(dx_in, R_in)          # inner boundary
# Extrude the ring region (start with 0; adjust if wrong)
sketch_and_extrude(face_pick_index=0, length_mm=L_extrude)

print("Done: conductors, cores, filler(%s), and shield created as separate bodies." % filler_mode)
