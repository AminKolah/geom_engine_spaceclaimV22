# Python Script, API Version = V252
import sys
from imp import reload

ROOT = r"C:\Users\kolahdouze\dev\geom_engine_spaceclaimV22"
if ROOT in sys.path:
    sys.path.remove(ROOT)
sys.path.insert(0, ROOT)

for m in [
    "engine",
    "engine.params",
    "engine.sketch_ops",
    "engine.body_ops",
    "engine.cable_build",
    "engine.naming",
    "engine.topology",
    "engine.surface_naming",
]:
    if m in sys.modules:
        del sys.modules[m]

from engine import params, sketch_ops, body_ops, cable_build, naming, topology, surface_naming

reload(params)
reload(sketch_ops)
reload(body_ops)
reload(cable_build)
reload(naming)
reload(topology)
reload(surface_naming)

from SpaceClaim.Api.V252 import *
from SpaceClaim.Api.V252.Geometry import *
from SpaceClaim.Api.V252.Modeler import *

ClearAll()

cfg = params.load_config(
    Parameters,
    input_file=r"C:\Users\kolahdouze\dev\geom_engine_spaceclaimV22\cases\input_case1.json",
    write_logs=True
)

api = {
    "GetRootPart": GetRootPart,
    "Delete": Delete,
    "BodySelection": BodySelection,
    "CurveSelection": CurveSelection,
    "FaceSelection": FaceSelection,
    "Fill": Fill,
    "Selection": Selection,
    "ExtrudeFaceOptions": ExtrudeFaceOptions,
    "ExtrudeType": ExtrudeType,
    "ExtrudeFaces": ExtrudeFaces,
    "Direction": Direction,
    "MM": MM,
    "ViewHelper": ViewHelper,
    "InteractionMode": InteractionMode,
    "Plane": Plane,
    "Frame": Frame,
    "Point": Point,
    "SketchCircle": SketchCircle,
    "SketchArc": SketchArc,
    "SketchLine": SketchLine,
    "SketchNurbs": SketchNurbs,
    "Point2D": Point2D,
    "Combine": Combine,
    "MakeSolidsOptions": MakeSolidsOptions,
    "NamedSelection": NamedSelection,
    "Window": Window,
    "Part": Part,
    "Component": Component,
    "ComponentHelper": ComponentHelper,
}

# --------------------------------
# 1) Build cable
# --------------------------------
bodies = cable_build.build_cable(api, cfg)

# --------------------------------
# 2) Move cable bodies into component
# --------------------------------
cable_body_names = [
    "conductor[1]",
    "conductor[2]",
    "single_core[1]",
    "single_core[2]",
    "single_core_merged",
    "Second_Extrusion",
    "Shield",
    "drain[1]",
    "drain[2]",
    "Overwrap",
]

cable_bodies_list = []
for n in cable_body_names:
    b = bodies.get(n)
    if b is not None:
        cable_bodies_list.append(b)

cable_comp = body_ops.get_or_create_component(
    GetRootPart, Window, Part, Component, "cable_bodies"
)

body_ops.move_bodies_to_component(
    ComponentHelper, BodySelection, cable_bodies_list, cable_comp
)

# Re-fetch after move
for n in cable_body_names:
    bodies[n] = body_ops.find_body_by_name_anywhere(GetRootPart, n)

# --------------------------------
# 3) Shared topology for bonded interfaces
# --------------------------------
topo_results = topology.apply_bonded_topology(
    api, cfg, bodies,
    verbose=True
)

# --------------------------------
# 4) Body-level named selections
# --------------------------------
ns_created = naming.create_all_body_named_selections(api, cfg, bodies)

# --------------------------------
# 5) Surface-level named selections
# --------------------------------
surface_ns_created = surface_naming.create_all_surface_named_selections(api, cfg, bodies)

print("Cable build complete.")

print("\n=== TOPOLOGY ===")
for k in sorted(topo_results.keys()):
    print("{:<25} {}".format(k, topo_results[k]))

print("\n=== BODY NAMED SELECTIONS ===")
for k in sorted(ns_created.keys()):
    print("{:<25} {}".format(k, ns_created[k]))

print("\n=== SURFACE NAMED SELECTIONS ===")
for k in sorted(surface_ns_created.keys()):
    print("{:<30} {}".format(k, surface_ns_created[k]))