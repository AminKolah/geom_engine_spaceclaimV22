# engine/section_props.py
from engine import sc_shims

def compute_stiffness_ratios(p, bodies):
    """
    Later: compute A, Ix, Iy (and EI_x/EI_y) from cross-section.
    For now: return placeholders.
    """
    return {"Ix": None, "Iy": None, "Ix_Iy": None}
