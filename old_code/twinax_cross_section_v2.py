# Python Script, API Version = V22

from SpaceClaim.Api.V22 import *
from SpaceClaim.Api.V22.Geometry import *

# -----------------------------
# INPUTS (mm)
# -----------------------------
D_cond    = 1.0
C2C       = 3.0

D_core    = C2C
D_core_eff = D_core
core_mode = "tangent"   # "tangent" or "explicit"
gap_core  = 0.0         # mm; used only in tangent mode

t_shield  = 0.10
W_outer   = 8.0
H_outer   = 5.0

# NEW: filler behavior
filler_mode = "fill"    # "fill" or "shell"

t_shell = W_outer/2 - D_core - t_shield
L_extrude = 200.0

# -----------------------------
# Setup plane + helpers
# -----------------------------
ViewHelper.SetSketchPlane(Plane.PlaneXY)

def P2(x, y):
    return Point2D.Create(MM(x), MM(y))

def draw_doubleD(dx, R):
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

    # Top/bottom connectors (this is what makes it a "fill")
    SketchLine.Create(P2(-dx, +R), P2(+dx, +R))
    SketchLine.Create(P2(-dx, -R), P2(+dx, -R))

def draw_doubleD_shell(dx, R, t_shell):
    # Arcs only (no connectors) => "shell" / side lobes only
    origin = P2(+dx, 0.0)
    start  = P2(+dx, +R)
    end    = P2(+dx, -R)
    SketchArc.CreateSweepArc(origin, start, end, True)

    origin = P2(-dx, 0.0)
    start  = P2(-dx, +R)
    end    = P2(-dx, -R)
    SketchArc.CreateSweepArc(origin, start, end, False)
    # Top/bottom connectors (this is what makes it a "fill")
    SketchLine.Create(P2(-dx, +R), P2(+dx, +R))
    SketchLine.Create(P2(-dx, -R), P2(+dx, -R))
    SketchLine.Create(P2(-dx, +R - t_shell), P2(+dx, +R - t_shell))
    SketchLine.Create(P2(-dx, -R + t_shell), P2(+dx, -R + t_shell))

def draw_doubleD_shell_consistent(W, H, t):
    # outer
    R  = H / 2.0
    dx = (W - H) / 2.0
    draw_doubleD(dx, R)

    # inner (same shape, offset inward)
    Wi = W - 2.0 * t
    Hi = H - 2.0 * t
    if Wi <= 0 or Hi <= 0:
        raise Exception("Shell thickness too large for given W/H.")
    Ri  = Hi / 2.0
    dxi = (Wi - Hi) / 2.0
    draw_doubleD(dxi, Ri)
    
# -----------------------------
# Derived
# -----------------------------
r_cond = D_cond / 2.0

if core_mode.lower() != "tangent":
    D_core_eff = C2C - gap_core/2

r_core = D_core_eff / 2.0

dx_cond = C2C / 2.0

# Outer shield
R_shield  = H_outer / 2.0
dx_shield = (W_outer - H_outer) / 2.0

# Inner filler (inside shield) inferred by offsetting outer by shield thickness
W_filler = W_outer - 2.0 * t_shield
H_filler = H_outer - 2.0 * t_shield
if W_filler <= 0 or H_filler <= 0:
    raise Exception("Shield thickness too large vs overall dimensions.")
if W_outer < H_outer:
    raise Exception("Invalid outer dimensions: W_outer must be >= H_outer for a double-D (two semicircles + straight).")

R_filler  = H_filler / 2.0
dx_filler = (W_filler - H_filler) / 2.0

# -----------------------------
# (A) Conductors + dielectric cores
# -----------------------------
SketchCircle.Create(P2(-dx_cond, 0.0), MM(r_cond))
SketchCircle.Create(P2(-dx_cond, 0.0), MM(r_core))

SketchCircle.Create(P2(+dx_cond, 0.0), MM(r_cond))
SketchCircle.Create(P2(+dx_cond, 0.0), MM(r_core))

# -----------------------------
# (B) Second extrusion: Fill vs Shell
# -----------------------------
mode = filler_mode.lower()
if mode == "fill":
    draw_doubleD(dx_shield, R_filler)
elif mode == "shell":
    draw_doubleD_shell(dx_filler, R_filler, t_shell)
else:
    raise Exception("Unknown filler_mode: %s (use 'fill' or 'shell')" % filler_mode)

# -----------------------------
# (C) Outer shield boundary
# -----------------------------
draw_doubleD(dx_cond, R_shield)

print("Cross-section sketch created. Filler mode =", filler_mode)
