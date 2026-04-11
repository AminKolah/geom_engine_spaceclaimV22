# -*- coding: utf-8 -*-

from SpaceClaim.Api.V252 import *
from SpaceClaim.Api.V252.Geometry import *
from SpaceClaim.Api.V252.Modeler import *

import math


# ------------------------------------------------------------
# basic utilities
# ------------------------------------------------------------

def _safe_name(obj):
    if obj is None:
        return "<None>"
    try:
        return str(obj.Name)
    except:
        pass
    try:
        return str(obj.GetName())
    except:
        pass
    return "<Unknown>"


def _delete_named_selection_if_exists(GetRootPart, Delete, Selection, ns_name):
    root = GetRootPart()
    try:
        for ns in list(root.NamedSelections):
            try:
                nm = ns.GetName()
            except:
                nm = getattr(ns, "Name", "")
            if nm == ns_name:
                try:
                    Delete.Execute(Selection.Create(ns))
                except:
                    try:
                        ns.Delete()
                    except:
                        pass
                return
    except:
        pass


def _create_face_ns(GetRootPart, Delete, Selection, FaceSelection, NamedSelection, faces, ns_name):
    if not faces:
        return None

    _delete_named_selection_if_exists(GetRootPart, Delete, Selection, ns_name)

    sel = FaceSelection.Create(faces)
    res = NamedSelection.Create(sel, Selection.Empty())

    try:
        res.CreatedNamedSelection.SetName(ns_name)
    except:
        try:
            res.CreatedNamedSelection.Name = ns_name
        except:
            return None

    return res.CreatedNamedSelection


def _face_repr_point(face):
    """
    Return a representative point guaranteed (as much as possible)
    to lie on the face.

    Strategy:
    1) take face bbox center
    2) project it onto face.Shape.Geometry
    3) fallback to first vertex position
    """
    try:
        bb = face.Shape.GetBoundingBox(Matrix.Identity)
        c = bb.Center

        # project bbox center onto the actual surface
        proj = face.Shape.Geometry.ProjectPoint(c)
        return proj.Point
    except:
        pass

    try:
        verts = list(face.Vertices)
        if verts:
            return verts[0].Position
    except:
        pass

    return None


def _body_bbox(body):
    return body.Shape.GetBoundingBox(Matrix.Identity)


def _is_endcap_by_z(face, body_bbox, z_tol_mm=0.01):
    """
    End caps are identified by face representative point lying very near
    the body's z-min or z-max.
    """
    p = _face_repr_point(face)
    if p is None:
        return False

    z_tol = z_tol_mm / 1000.0
    z = float(p.Z)
    cz = float(body_bbox.Center.Z)
    hz = 0.5 * float(body_bbox.Size.Z)

    zmin = cz - hz
    zmax = cz + hz

    return (abs(z - zmin) < z_tol) or (abs(z - zmax) < z_tol)


def _rho_local(x, y, cx, cy):
    dx = x - cx
    dy = y - cy
    return math.sqrt(dx * dx + dy * dy)


def _doubleD_boundary_radius(x, y, cx, dx_half, R):
    """
    For a Double-D centered around y=0 with left/right arc centers at:
      (cx - dx_half, 0), (cx + dx_half, 0)

    return a measure of distance to the nearest circular/flat boundary.
    This is an approximate classifier, not an exact signed distance.
    """
    xr = x - cx

    # top/bottom flats region
    if abs(xr) <= dx_half:
        # closeness to flat is based on |y|
        return abs(y)

    # circular caps
    if xr > dx_half:
        ccx = cx + dx_half
    else:
        ccx = cx - dx_half

    return math.sqrt((x - ccx) ** 2 + y ** 2)


def _elliptic_like_radius(x, y, cx, cy):
    """
    Simple radial proxy from a center; used only as fallback
    if you decide to support elliptic later in this classifier.
    """
    return math.sqrt((x - cx) ** 2 + (y - cy) ** 2)


# ------------------------------------------------------------
# conductor classification
# ------------------------------------------------------------

def _classify_conductor_faces(body, center_x, center_y, r_cond):
    """
    For a conductor:
      - remove end caps by z-extrema
      - everything else is outer
    """
    bbox = _body_bbox(body)
    outer_faces = []

    for f in list(body.Faces):
        if _is_endcap_by_z(f, bbox):
            continue
        outer_faces.append(f)

    return outer_faces


# ------------------------------------------------------------
# single-core classification
# ------------------------------------------------------------

def _classify_single_core_faces(body, center_x, center_y, r_cond, r_core,
                                touching_x_sign,
                                radial_tol_mm=0.02,
                                touching_y_tol_mm=0.02,
                                touching_xcos=0.95):
    """
    Buckets:
      - outer_faces: external insulation surface, near local radius r_core
      - conductor_faces: conductor-hole inner surface, near local radius r_cond
      - lumen_faces: other interior cavity walls
      - touching_faces: tiny interface face between the two touching cores

    touching_x_sign:
      +1 means touching face normal-ish direction expected toward +X
      -1 means toward -X
    """
    bbox = _body_bbox(body)

    outer_faces = []
    conductor_faces = []
    lumen_faces = []
    touching_faces = []

    radial_tol = radial_tol_mm / 1000.0
    touching_y_tol = touching_y_tol_mm / 1000.0

    for f in list(body.Faces):
        if _is_endcap_by_z(f, bbox):
            continue

        p = _face_repr_point(f)
        if p is None:
            continue

        x = float(p.X)
        y = float(p.Y)

        # touching face heuristic:
        # near global midline y ~ 0, and x close to the inter-core interface region
        # We still allow a mild use of planar normal if available.
        is_touching = False
        try:
            n = f.Plane.Normal
            nx = float(n.X)
            ny = float(n.Y)
            if abs(nx) > touching_xcos and abs(ny) < 0.2 and abs(y - center_y) < touching_y_tol:
                if touching_x_sign > 0 and nx > 0:
                    is_touching = True
                elif touching_x_sign < 0 and nx < 0:
                    is_touching = True
        except:
            pass

        if is_touching:
            touching_faces.append(f)
            continue

        rho = _rho_local(x, y, center_x, center_y)

        if abs(rho - r_core / 1000.0) < radial_tol:
            outer_faces.append(f)
        elif abs(rho - r_cond / 1000.0) < radial_tol:
            conductor_faces.append(f)
        else:
            lumen_faces.append(f)

    return outer_faces, conductor_faces, lumen_faces, touching_faces


# ------------------------------------------------------------
# layer classification (Second_Extrusion, Shield, Overwrap)
# ------------------------------------------------------------

def _classify_layer_faces_doubleD(body, dx_outer, R_outer, dx_inner, R_inner,
                                  center_x=0.0, center_y=0.0,
                                  radial_tol_mm=0.03):
    """
    Construction-aware classification for Double-D layers.

    For each side face:
      - compute representative point
      - compare closeness to the outer boundary vs inner boundary
      - whichever boundary it is closer to determines outer/inner

    This avoids normal guessing entirely.
    """
    bbox = _body_bbox(body)

    outer_faces = []
    inner_faces = []

    tol = radial_tol_mm / 1000.0

    for f in list(body.Faces):
        if _is_endcap_by_z(f, bbox):
            continue

        p = _face_repr_point(f)
        if p is None:
            continue

        x = float(p.X)
        y = float(p.Y)

        ro = _doubleD_boundary_radius(x, y, center_x, dx_outer / 1000.0, R_outer / 1000.0)
        ri = _doubleD_boundary_radius(x, y, center_x, dx_inner / 1000.0, R_inner / 1000.0)

        # Face closer to outer or inner reference boundary?
        # We compare the mismatch to the nominal boundary size.
        # Smaller mismatch wins.
        do = abs(ro - (R_outer / 1000.0 if abs(x - center_x) > dx_outer / 1000.0 else abs(y)))
        di = abs(ri - (R_inner / 1000.0 if abs(x - center_x) > dx_inner / 1000.0 else abs(y)))

        if do + tol < di:
            outer_faces.append(f)
        elif di + tol < do:
            inner_faces.append(f)
        else:
            # fallback: farther from center usually means outer
            rho = _rho_local(x, y, center_x, center_y)
            # use rough proxy based on average of inner/outer envelopes
            if rho > 0.5 * ((R_outer / 1000.0) + (R_inner / 1000.0)):
                outer_faces.append(f)
            else:
                inner_faces.append(f)

    return outer_faces, inner_faces


# ------------------------------------------------------------
# public creator
# ------------------------------------------------------------

def create_all_surface_named_selections(api, cfg, bodies):
    GetRootPart = api["GetRootPart"]
    Delete = api["Delete"]
    Selection = api["Selection"]
    FaceSelection = api["FaceSelection"]
    NamedSelection = api["NamedSelection"]

    created = {}

    # --------------------------
    # Conductors
    # --------------------------
    c1 = bodies.get("conductor[1]")
    if c1 is not None:
        faces = _classify_conductor_faces(
            c1,
            center_x=-cfg["dx_cond"] / 1000.0,
            center_y=0.0,
            r_cond=cfg["D_cond"]
        )
        created["NS_conductor[1]_Outer"] = _create_face_ns(
            GetRootPart, Delete, Selection, FaceSelection, NamedSelection,
            faces, "NS_conductor[1]_Outer"
        )
    else:
        created["NS_conductor[1]_Outer"] = None

    c2 = bodies.get("conductor[2]")
    if c2 is not None:
        faces = _classify_conductor_faces(
            c2,
            center_x=cfg["dx_cond"] / 1000.0,
            center_y=0.0,
            r_cond=cfg["D_cond"]
        )
        created["NS_conductor[2]_Outer"] = _create_face_ns(
            GetRootPart, Delete, Selection, FaceSelection, NamedSelection,
            faces, "NS_conductor[2]_Outer"
        )
    else:
        created["NS_conductor[2]_Outer"] = None

    # --------------------------
    # Single cores
    # --------------------------
    if cfg["core_mode"] == "separate":
        s1 = bodies.get("single_core[1]")
        if s1 is not None:
            outer_faces, conductor_faces, lumen_faces, touching_faces = _classify_single_core_faces(
                s1,
                center_x=-cfg["dx_cond"] / 1000.0,
                center_y=0.0,
                r_cond=0.5 * cfg["D_cond"],
                r_core=0.5 * cfg["D_core"],
                touching_x_sign=+1
            )
            created["NS_single_core[1]_Outer"] = _create_face_ns(
                GetRootPart, Delete, Selection, FaceSelection, NamedSelection,
                outer_faces, "NS_single_core[1]_Outer"
            )
            created["NS_single_core[1]_Conductor"] = _create_face_ns(
                GetRootPart, Delete, Selection, FaceSelection, NamedSelection,
                conductor_faces, "NS_single_core[1]_Conductor"
            )
            created["NS_single_core[1]_Lumens"] = _create_face_ns(
                GetRootPart, Delete, Selection, FaceSelection, NamedSelection,
                lumen_faces, "NS_single_core[1]_Lumens"
            )
            created["NS_single_core[1]_Touching"] = _create_face_ns(
                GetRootPart, Delete, Selection, FaceSelection, NamedSelection,
                touching_faces, "NS_single_core[1]_Touching"
            )
        else:
            created["NS_single_core[1]_Outer"] = None
            created["NS_single_core[1]_Conductor"] = None
            created["NS_single_core[1]_Lumens"] = None
            created["NS_single_core[1]_Touching"] = None

        s2 = bodies.get("single_core[2]")
        if s2 is not None:
            outer_faces, conductor_faces, lumen_faces, touching_faces = _classify_single_core_faces(
                s2,
                center_x=cfg["dx_cond"] / 1000.0,
                center_y=0.0,
                r_cond=cfg["D_cond"],
                r_core=cfg["D_core"],
                touching_x_sign=-1
            )
            created["NS_single_core[2]_Outer"] = _create_face_ns(
                GetRootPart, Delete, Selection, FaceSelection, NamedSelection,
                outer_faces, "NS_single_core[2]_Outer"
            )
            created["NS_single_core[2]_Conductor"] = _create_face_ns(
                GetRootPart, Delete, Selection, FaceSelection, NamedSelection,
                conductor_faces, "NS_single_core[2]_Conductor"
            )
            created["NS_single_core[2]_Lumens"] = _create_face_ns(
                GetRootPart, Delete, Selection, FaceSelection, NamedSelection,
                lumen_faces, "NS_single_core[2]_Lumens"
            )
            created["NS_single_core[2]_Touching"] = _create_face_ns(
                GetRootPart, Delete, Selection, FaceSelection, NamedSelection,
                touching_faces, "NS_single_core[2]_Touching"
            )
        else:
            created["NS_single_core[2]_Outer"] = None
            created["NS_single_core[2]_Conductor"] = None
            created["NS_single_core[2]_Lumens"] = None
            created["NS_single_core[2]_Touching"] = None

    else:
        created["NS_single_core_merged_Outer"] = None
        created["NS_single_core_merged_Conductor"] = None
        created["NS_single_core_merged_Lumens"] = None
        created["NS_single_core_merged_Touching"] = None

    # --------------------------
    # Second_Extrusion
    # --------------------------
    second_extrusion = bodies.get("Second_Extrusion")
    if second_extrusion is not None:
        if cfg["is_elliptic"]:
            # fallback for now
            outer_faces = []
            inner_faces = []
        else:
            outer_faces, inner_faces = _classify_layer_faces_doubleD(
                second_extrusion,
                dx_outer=cfg["dx_in"],
                R_outer=cfg["R_in"],
                dx_inner=cfg["dx_cond"],
                R_inner=cfg["D_core"] / 2.0
            )
        created["NS_Second_Extrusion_Outer"] = _create_face_ns(
            GetRootPart, Delete, Selection, FaceSelection, NamedSelection,
            outer_faces, "NS_Second_Extrusion_Outer"
        )
        created["NS_Second_Extrusion_Inner"] = _create_face_ns(
            GetRootPart, Delete, Selection, FaceSelection, NamedSelection,
            inner_faces, "NS_Second_Extrusion_Inner"
        )
    else:
        created["NS_Second_Extrusion_Outer"] = None
        created["NS_Second_Extrusion_Inner"] = None

    # --------------------------
    # Shield
    # --------------------------
    shield = bodies.get("Shield")
    if shield is not None:
        if cfg["is_elliptic"]:
            outer_faces = []
            inner_faces = []
        else:
            outer_faces, inner_faces = _classify_layer_faces_doubleD(
                shield,
                dx_outer=cfg["dx_shield"],
                R_outer=cfg["R_shield"],
                dx_inner=cfg["dx_in"],
                R_inner=cfg["R_in"]
            )
        created["NS_Shield_Outer"] = _create_face_ns(
            GetRootPart, Delete, Selection, FaceSelection, NamedSelection,
            outer_faces, "NS_Shield_Outer"
        )
        created["NS_Shield_Inner"] = _create_face_ns(
            GetRootPart, Delete, Selection, FaceSelection, NamedSelection,
            inner_faces, "NS_Shield_Inner"
        )
    else:
        created["NS_Shield_Outer"] = None
        created["NS_Shield_Inner"] = None

    # --------------------------
    # Overwrap
    # --------------------------
    overwrap = bodies.get("Overwrap")
    if overwrap is not None:
        if cfg["is_elliptic"]:
            outer_faces = []
            inner_faces = []
        else:
            outer_faces, inner_faces = _classify_layer_faces_doubleD(
                overwrap,
                dx_outer=cfg["dx_shield"],
                R_outer=cfg["R_shield"] + cfg["t_overwrap"],
                dx_inner=cfg["dx_shield"],
                R_inner=cfg["R_shield"]
            )
        created["NS_Overwrap_Outer"] = _create_face_ns(
            GetRootPart, Delete, Selection, FaceSelection, NamedSelection,
            outer_faces, "NS_Overwrap_Outer"
        )
        created["NS_Overwrap_Inner"] = _create_face_ns(
            GetRootPart, Delete, Selection, FaceSelection, NamedSelection,
            inner_faces, "NS_Overwrap_Inner"
        )
    else:
        created["NS_Overwrap_Outer"] = None
        created["NS_Overwrap_Inner"] = None

    return created