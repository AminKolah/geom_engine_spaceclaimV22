# -*- coding: utf-8 -*-

from SpaceClaim.Api.V252 import *
from SpaceClaim.Api.V252.Geometry import *
from SpaceClaim.Api.V252.Modeler import *


def _safe_name(obj):
    if obj is None:
        return "<None>"
    try:
        n = getattr(obj, "Name", None)
        if n is not None:
            return str(n)
    except:
        pass
    try:
        if hasattr(obj, "GetName"):
            return str(obj.GetName())
    except:
        pass
    return "<Unknown>"


def _get_parent(body):
    try:
        return body.Parent
    except:
        return None


def _set_shared_topology_on_parent(body, verbose=True):
    if body is None:
        if verbose:
            print("SharedTopology skipped: body is None")
        return False

    parent = _get_parent(body)
    if parent is None:
        if verbose:
            print("SharedTopology skipped: no parent for", _safe_name(body))
        return False

    try:
        if hasattr(parent, "SharedTopology"):
            parent.SharedTopology = SharedTopology.Share
            if verbose:
                print("SharedTopology set to Share on parent of", _safe_name(body))
            return True
    except Exception as e:
        if verbose:
            print("Failed setting SharedTopology on parent of", _safe_name(body), ":", e)
        return False

    if verbose:
        print("Parent has no SharedTopology attribute for", _safe_name(body))
    return False


def share_topology_pair(body_a, body_b, verbose=True):
    """
    Minimal version: only sets SharedTopology = Share on the parents
    of the two bodies. No imprint/repair operations.
    """
    if body_a is None or body_b is None:
        if verbose:
            print("share_topology_pair skipped:", _safe_name(body_a), _safe_name(body_b))
        return False

    ok_a = _set_shared_topology_on_parent(body_a, verbose=verbose)
    ok_b = _set_shared_topology_on_parent(body_b, verbose=verbose)

    if verbose:
        print("Requested shared topology between", _safe_name(body_a), "and", _safe_name(body_b))

    return (ok_a or ok_b)


def apply_bonded_topology(api, cfg, bodies, verbose=True):
    """
    Applies shared topology only to the interfaces intended to be bonded/conformal:
      - conductor[1] <-> single_core[1]
      - conductor[2] <-> single_core[2]
      - single_core(s) <-> Second_Extrusion

    This version only uses SharedTopology on parent/components.
    """
    results = {}

    core_mode = cfg["core_mode"]

    c1 = bodies.get("conductor[1]")
    c2 = bodies.get("conductor[2]")
    filler = bodies.get("Second_Extrusion")

    if core_mode == "separate":
        s1 = bodies.get("single_core[1]")
        s2 = bodies.get("single_core[2]")

        results["conductor1_core1"] = share_topology_pair(c1, s1, verbose=verbose)
        results["conductor2_core2"] = share_topology_pair(c2, s2, verbose=verbose)
        results["core1_filler"] = share_topology_pair(s1, filler, verbose=verbose)
        results["core2_filler"] = share_topology_pair(s2, filler, verbose=verbose)

    else:
        s_merged = bodies.get("single_core_merged")

        results["conductor1_coreMerged"] = share_topology_pair(c1, s_merged, verbose=verbose)
        results["conductor2_coreMerged"] = share_topology_pair(c2, s_merged, verbose=verbose)
        results["coreMerged_filler"] = share_topology_pair(s_merged, filler, verbose=verbose)

    return results