# -*- coding: utf-8 -*-



import math

from SpaceClaim.Api.V252 import *               # <-- IMPORTANT (ViewHelper, InteractionMode, etc.)
from SpaceClaim.Api.V252.Geometry import *
from SpaceClaim.Api.V252.Modeler import *
from engine import sc_shims
# -----------------------------------------------------------------------------
# View / sketch solidify
# -----------------------------------------------------------------------------

def solidify_sketch():
    sc_shims.set_view_mode_solid()

# -----------------------------------------------------------------------------
# Face selection helpers
# -----------------------------------------------------------------------------

def _largest_xy_face(body):
    """
    Return the largest planar-ish face (by Area) on a body.
    Used by extrude helpers that need a start face.
    """
    faces = list(body.Faces)
    if not faces:
        return None
    return max(faces, key=lambda f: f.Area)

# -----------------------------------------------------------------------------
# Extrude helpers (lifted from monolith)
# -----------------------------------------------------------------------------

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

def union_bodies(target, tool):
    """
    Boolean union tool into target. Returns unioned target body.
    """
    sel = Selection.Create(target)
    opts = CombineOptions()
    opts.Operation = CombineOperation.Union
    Combine.Execute(sel, Selection.Create(tool), opts)
    return target

def cleanup_sheet_bodies(component=None):
    """
    Delete sheet bodies (surface bodies) under root or a component.
    """
    part = GetRootPart()
    bodies = []
    if component is None:
        bodies = list(part.Bodies)
    else:
        bodies = list(component.GetAllBodies())

    for b in bodies:
        try:
            if b.Shape and b.Shape.ShapeType == ShapeType.SheetBody:
                b.Delete()
        except:
            pass

# -----------------------------------------------------------------------------
# Component helpers (monolith utilities)
# -----------------------------------------------------------------------------

def _comp_by_name(name):
    for c in GetRootPart().Components:
        if c.Name == name:
            return c
    return None

def get_or_create_component(name):
    c = _comp_by_name(name)
    if c is not None:
        return c
    c = Component.Create(GetRootPart(), name).CreatedComponent
    c.Name = name
    return c

def move_body_to_component(body, component):
    """
    Move an existing body into a component.
    """
    if body is None or component is None:
        return
    try:
        ComponentHelper.MoveBodiesToComponent(Selection.Create(body), component)
    except:
        # older API variants might require different calls; we'll patch later if needed
        ComponentHelper.MoveBodiesToComponent(Selection.Create(body), component)

def move_body_to_component_once(body, component):
    """
    Avoid double-moving if body already lives under that component.
    """
    try:
        if body.Parent == component:
            return
    except:
        pass
    move_body_to_component(body, component)

def move_all_root_bodies_to_component(component_name):
    """
    Collect all root bodies into a component (common in your monolith).
    """
    comp = get_or_create_component(component_name)
    root = GetRootPart()
    for b in list(root.Bodies):
        try:
            move_body_to_component_once(b, comp)
        except:
            pass
    return comp

# -----------------------------------------------------------------------------
# Volume utilities (used in your monolith cleanup)
# -----------------------------------------------------------------------------

def body_volume_mm3(body):
    try:
        mp = MeasureHelper.GetMassProperties(Selection.Create(body))
        return mp.Volume
    except:
        return 0.0

def delete_zero_volume_bodies_in_component(component_name, eps_vol=1e-9):
    comp = _comp_by_name(component_name)
    if comp is None:
        return
    for b in list(comp.GetAllBodies()):
        try:
            if body_volume_mm3(b) <= eps_vol:
                b.Delete()
        except:
            pass

def delete_bodies_by_exact_name_in_component(component_name, body_name):
    comp = _comp_by_name(component_name)
    if comp is None:
        return
    for b in list(comp.GetAllBodies()):
        try:
            if b.Name == body_name:
                b.Delete()
        except:
            pass

# -----------------------------------------------------------------------------
# Shared topology helpers (ported from monolith; keep BOTH working + WIP)
# -----------------------------------------------------------------------------

def _safe_name(obj):
    try:
        return obj.Name
    except:
        try:
            return str(obj)
        except:
            return "<unnamed>"

def _set_shared_topology_on_parent(body, verbose=True):
    """
    Sets SharedTopology on the parent component if the property exists.
    Returns True if set succeeded, else False.
    """
    if body is None:
        return False
    try:
        comp = body.Parent
        if hasattr(comp, 'SharedTopology'):
            comp.SharedTopology = SharedTopology.Share
            if verbose:
                print("Component Property set to 'Share' for:", _safe_name(body))
            return True
    except Exception as e:
        if verbose:
            print("Failed to set SharedTopology property for", _safe_name(body), ":", e)
    return False


# -------------------------
# WORKING VERSION (keep as default)
# -------------------------
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

def share_topology_group(bodies, do_imprint=False, imprint_tol_mm=None, verbose=True):
    bodies = [b for b in bodies if b is not None]
    if len(bodies) < 2:
        return

    try:
        for b in bodies:
            comp = getattr(b, "Parent", None)
            if comp is not None and hasattr(comp, "SharedTopology"):
                comp.SharedTopology = SharedTopology.Share
        if verbose:
            print("SharedTopology set to Share on parent(s) (best effort).")
    except Exception as e:
        if verbose:
            print("Failed setting SharedTopology:", e)

    if not do_imprint:
        return

    selection = BodySelection.Create(bodies)
    try:
        opts = ImprintOptions()
        if imprint_tol_mm is not None and hasattr(opts, "Tolerance"):
            opts.Tolerance = MM(imprint_tol_mm)
        Imprint.Execute(selection, opts)
        if verbose:
            print("Imprint successful on group of", len(bodies), "bodies.")
        return
    except Exception as e_imprint:
        if verbose:
            print("Imprint failed:", e_imprint)

    try:
        Share.Share(selection)
        if verbose:
            print("Repair.Share successful on group.")
    except Exception as e_share:
        if verbose:
            print("Repair.Share failed:", e_share)

# -------------------------
# WIP VERSION (keep for later troubleshooting)
# -------------------------
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
            print("share_topology_pair_WIP: same body provided; skipping:", _safe_name(body_a))
        return True

    # 1) Always do SharedTopology flag first (low-risk, Mechanical-friendly)
    ok_a = _set_shared_topology_on_parent(body_a, verbose=verbose)
    ok_b = _set_shared_topology_on_parent(body_b, verbose=verbose)

    # If user only wants the “WB share topology” behavior, stop here.
    if mode.lower().strip() in ["share_only", "share"]:
        if verbose:
            print("share_topology_pair_WIP: mode=share_only; done.")
        return (ok_a or ok_b)

    # 2) Optional: physical imprint / repair-share (higher risk for HFSS)
    sel = None
    try:
        sel = BodySelection.Create([body_a, body_b])
    except Exception as e:
        if verbose:
            print("share_topology_pair_WIP: failed to create BodySelection:", e)
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
        print("share_topology_pair_WIP: all methods failed for:",
              _safe_name(body_a), "<->", _safe_name(body_b))
    return False
