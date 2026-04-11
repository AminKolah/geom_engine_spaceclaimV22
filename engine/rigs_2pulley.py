# engine/rigs_2pulley.py
from engine import sc_shims

def build(p, bodies):
    """
    Build a two-pulley rig around the cable geometry.

    Inputs:
      p: Params
      bodies: dict from cable_build.build(p)
    Outputs:
      Optionally returns a dict of rig bodies.
    """
    sc_shims.log("[rigs_2pulley] build() start")

    rig = {}

    cable = bodies.get("Cable", None) or bodies.get("Second_Extrusion", None)
    if cable is None:
        sc_shims.log("[rigs_2pulley] No cable-like body found in 'bodies'. Skipping rig creation.")
        return rig

    pulley_radius = getattr(p, "pulley_radius", 10.0)
    pulley_width  = getattr(p, "pulley_width",  5.0)
    pulley_c2c    = getattr(p, "pulley_c2c",   40.0)

    groove_clear  = getattr(p, "groove_clear", 0.2)

    sc_shims.log("[rigs_2pulley] R={} W={} C2C={} clear={}".format(
        pulley_radius, pulley_width, pulley_c2c, groove_clear
    ))

    # TODO: rig["Pulley_1"] = create_pulley(...)
    # TODO: rig["Pulley_2"] = create_pulley(...)
    # TODO: rig["GrooveTool"] = create_groove_tool(...); subtract...

    sc_shims.log("[rigs_2pulley] build() done")
    return rig


# -----------------------------------------------------------------------------
# Optional helper function placeholders (fill in with your existing code)
# -----------------------------------------------------------------------------

def create_pulley(name, radius, width, center_point, axis_dir):
    """
    Create a pulley wheel (likely via sketch + revolve).
    """
    raise NotImplementedError

def create_groove_tool(name, width, cable_W, cable_H, clear, center_y, center_z):
    """
    Create groove cutter tool (your 'extrude along x' version can live here).
    """
    raise NotImplementedError
