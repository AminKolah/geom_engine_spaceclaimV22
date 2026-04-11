# -*- coding: utf-8 -*-
# engine/sc_shims.py

import sys

def _get_host_module():
    """
    SpaceClaim Script Editor can run under different host module names.
    We try several likely names.
    """
    for name in ("SpaceClaim_Script", "__main__", "SpaceClaimScript", "SpaceClaim"):
        m = sys.modules.get(name, None)
        if m is not None:
            return m
    return None

_host = _get_host_module()

# Best-effort API imports (may or may not expose UI helpers like ViewHelper)
try:
    from SpaceClaim.Api.V252 import *
    from SpaceClaim.Api.V252.Geometry import *
    from SpaceClaim.Api.V252.Modeler import *
except:
    pass

def _bind(name):
    """
    Prefer host-injected symbol; fall back to whatever imports provided.
    """
    try:
        if _host is not None and hasattr(_host, name):
            return getattr(_host, name)
    except:
        pass

    if name in globals():
        return globals()[name]

    return None

def get(name, default=None):
    v = _bind(name)
    return v if v is not None else default

# ---- Core shims ----
ViewHelper      = get("ViewHelper")
InteractionMode = get("InteractionMode")

Plane      = get("Plane")
Frame      = get("Frame")
Point      = get("Point")
Direction  = get("Direction")
MM         = get("MM")
Point2D    = get("Point2D")

GetRootPart    = get("GetRootPart")
Selection      = get("Selection")
FaceSelection  = get("FaceSelection")
BodySelection  = get("BodySelection")
Delete         = get("Delete")

# ---- Optional helpers (don’t fail if missing) ----
ViewModeHelper = get("ViewModeHelper")
SketchHelper   = get("SketchHelper")

# ---- Safe wrappers ----
def log(msg):
    try:
        print(msg)
    except:
        pass

def set_view_mode_solid():
    """
    Make view solid without hard dependency on ViewHelper.
    """
    try:
        if ViewHelper is not None and InteractionMode is not None:
            ViewHelper.SetViewMode(InteractionMode.Solid)
            return True
    except:
        pass

    try:
        if ViewModeHelper is not None and InteractionMode is not None:
            ViewModeHelper.SetViewMode(InteractionMode.Solid)
            return True
    except:
        pass

    return False

def set_sketch_plane(plane):
    """
    Set sketch plane without hard dependency on ViewHelper.
    """
    try:
        if ViewHelper is not None:
            ViewHelper.SetSketchPlane(plane)
            return True
    except:
        pass

    try:
        if SketchHelper is not None:
            SketchHelper.SetSketchPlane(plane)
            return True
    except:
        pass

    return False

# ---- Sketch primitives ----
def _first_attr(names):
    """Return the first non-None attribute found across known namespaces."""
    # 1) already in this module globals (from import *)
    for n in names:
        v = get(n)
        if v is not None:
            return v

    # 2) try common namespaces explicitly
    candidates = []
    try:
        import SpaceClaim.Api.V252 as _A
        candidates.append(_A)
    except Exception as e:
        log("sc_shims: cannot import SpaceClaim.Api.V252: %s" % e)

    try:
        import SpaceClaim.Api.V252.Modeler as _M
        candidates.append(_M)
    except Exception as e:
        log("sc_shims: cannot import SpaceClaim.Api.V252.Modeler: %s" % e)

    try:
        import SpaceClaim.Api.V252.Sketch as _S
        candidates.append(_S)
    except:
        # not all builds have this
        pass

    for mod in candidates:
        for n in names:
            try:
                if hasattr(mod, n):
                    return getattr(mod, n)
            except:
                pass

    # 3) last resort: scan loaded modules for the symbol
    for modname, mod in sys.modules.items():
        try:
            for n in names:
                if hasattr(mod, n):
                    return getattr(mod, n)
        except:
            pass

    return None

SketchCircle = _first_attr(["SketchCircle"])
SketchLine   = _first_attr(["SketchLine"])
SketchArc    = _first_attr(["SketchArc"])
SketchNurbs  = _first_attr(["SketchNurbs"])

log("sc_shims: SketchCircle=%s" % (SketchCircle,))
log("sc_shims: SketchLine=%s" % (SketchLine,))
log("sc_shims: SketchArc=%s" % (SketchArc,))
log("sc_shims: SketchNurbs=%s" % (SketchNurbs,))
