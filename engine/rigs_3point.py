# engine/rigs_3point.py
from engine import sc_shims

def build(p, bodies):
    """
    Build a 3-point bending rig around the cable geometry.

    Inputs:
      p: Params
      bodies: dict from cable_build.build(p)
              e.g. bodies["Cable"], bodies["Shield"], etc. (whatever you decide)
    Outputs:
      Optionally returns a dict of rig bodies.
    """
    sc_shims.log("[rigs_3point] build() start")

    rig = {}

    # -------------------------------------------------------------------------
    # TODO: move your 3-point bending rig code here
    #
    # Typical structure:
    #   1) locate cable body or reference envelope
    #   2) create supports (two cylinders or blocks) at +/- span/2
    #   3) create loading nose (cylinder) at midspan
    #   4) split/trim supports/nose for "good vs bad bending direction"
    #   5) name/organize components and named selections
    # -------------------------------------------------------------------------

    # Example placeholders (safe; won't crash if bodies are missing):
    cable = bodies.get("Cable", None) or bodies.get("Second_Extrusion", None)
    if cable is None:
        sc_shims.log("[rigs_3point] No cable-like body found in 'bodies'. Skipping rig creation.")
        return rig

    # You can keep these params in params.py eventually
    span = getattr(p, "span", None)
    if span is None:
        span = 50.0  # units: your model units

    nose_radius = getattr(p, "nose_radius", 2.0)
    support_radius = getattr(p, "support_radius", 2.0)

    sc_shims.log("[rigs_3point] span={} nose_R={} support_R={}".format(span, nose_radius, support_radius))

    # TODO: rig["Support_1"] = create_support(...)
    # TODO: rig["Support_2"] = create_support(...)
    # TODO: rig["Nose"]      = create_loading_nose(...)

    sc_shims.log("[rigs_3point] build() done")
    return rig


# -----------------------------------------------------------------------------
# Optional helper function placeholders (fill in with your existing code)
# -----------------------------------------------------------------------------

def create_support(name, radius, length, center_point, axis_dir):
    """
    Create a cylindrical support body.
    Implement using your existing sketch/extrude or primitive creation.
    """
    raise NotImplementedError

def create_loading_nose(name, radius, length, center_point, axis_dir):
    """
    Create a cylindrical loading nose body.
    """
    raise NotImplementedError
