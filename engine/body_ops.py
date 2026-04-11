


def find_body_by_name(GetRootPart, name):
    root = GetRootPart()
    for b in list(root.Bodies):
        if getattr(b, "Name", "") == name:
            return b
    for comp in list(root.Components):
        try:
            for b in list(comp.Content.Bodies):
                if getattr(b, "Name", "") == name:
                    return b
        except:
            pass
    return None


def _largest_xy_face(body):
    best = None
    bestA = -1.0
    for f in list(body.Faces):
        try:
            n = f.Plane.Normal
            if abs(n.Z) > 0.95:
                if f.Area > bestA:
                    best = f
                    bestA = f.Area
        except:
            pass
    return best


def extrude_and_name(GetRootPart, Delete, BodySelection, CurveSelection,
                     Fill, Selection, ExtrudeFaceOptions, ExtrudeType,ExtrudeFaces,
                     FaceSelection, Direction, MM, solidify_sketch,
                     clear_all_sketch_curves,
                     name, length_mm, pick_largest=True):
    root = GetRootPart()
    before = list(root.Bodies)

    solidify_sketch()

    after = list(root.Bodies)
    new_bodies = [b for b in after if b not in before]

    if not new_bodies:
        try:
            crvs = list(root.Curves)
            if not crvs:
                raise Exception("No sketch curves found to Fill for '%s'." % name)
            Fill.Execute(CurveSelection.Create(crvs))
        except Exception as e:
            raise Exception("Solidify produced NO sketch-fill body for '%s' and Fill failed: %s" % (name, e))

        after2 = list(root.Bodies)
        new_bodies = [b for b in after2 if b not in before]
        if not new_bodies:
            raise Exception("Solidify+Fill produced NO sketch-fill body for '%s'." % name)

    def best_xy_area(b):
        f = _largest_xy_face(b)
        if f is None:
            try:
                f = max(list(b.Faces), key=lambda ff: ff.Area)
            except:
                return -1.0
        try:
            return float(f.Area)
        except:
            return -1.0

    if pick_largest:
        temp_body = max(new_bodies, key=best_xy_area)
    else:
        bodies_sorted = sorted(new_bodies, key=best_xy_area)
        temp_body = bodies_sorted[len(bodies_sorted) // 2]

    face = _largest_xy_face(temp_body)
    if face is None:
        face = max(list(temp_body.Faces), key=lambda f: f.Area)

    options = ExtrudeFaceOptions()
    options.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirZ, MM(length_mm), options)

    created = list(res.CreatedBodies)
    if not created:
        raise Exception("Extrude produced no result for '%s'." % name)

    new_body = created[0]
    try:
        new_body.SetName(name)
    except:
        new_body.Name = name

    try:
        Delete.Execute(BodySelection.Create(new_bodies))
    except:
        pass

    try:
        clear_all_sketch_curves()
    except:
        pass

    return new_body


def extrude_from_explicit_curves(GetRootPart, Delete, BodySelection, CurveSelection,
                                 Fill, ExtrudeFaceOptions, ExtrudeType, ExtrudeFaces,
                                 FaceSelection, Direction, MM, ViewHelper,
                                 InteractionMode, clear_all_sketch_curves,
                                 name, length_mm, curves):
    root = GetRootPart()
    before = list(root.Bodies)

    if not curves:
        raise Exception("No curves provided for '%s'." % name)

    Fill.Execute(CurveSelection.Create(curves))

    after = list(root.Bodies)
    new_bodies = [b for b in after if b not in before]

    if not new_bodies:
        raise Exception("Explicit Fill created no body for '%s'." % name)

    body = new_bodies[-1]
    faces = list(body.Faces)
    if not faces:
        raise Exception("Filled body for '%s' has no faces." % name)

    face = max(faces, key=lambda f: f.Area)

    opts = ExtrudeFaceOptions()
    opts.ExtrudeType = ExtrudeType.ForceIndependent
    res = ExtrudeFaces.Execute(FaceSelection.Create(face), Direction.DirZ, MM(length_mm), opts)

    created = list(res.CreatedBodies)
    if not created:
        raise Exception("Extrude produced no result for '%s'." % name)

    new_body = created[0]
    try:
        new_body.SetName(name)
    except:
        new_body.Name = name

    try:
        Delete.Execute(BodySelection.Create(new_bodies))
    except:
        pass

    try:
        clear_all_sketch_curves()
    except:
        pass

    try:
        ViewHelper.SetViewMode(InteractionMode.Solid)
    except:
        pass

    return new_body


def union_bodies(Combine, Selection, name, bodies):
    bodies = [b for b in bodies if b is not None]
    if len(bodies) < 2:
        return bodies[0] if bodies else None

    try:
        res = Combine.Merge(Selection.Create(bodies))
        try:
            out = res.MergedBody
        except:
            try:
                out = list(res.CreatedBodies)[0]
            except:
                out = bodies[0]
        try:
            out.SetName(name)
        except:
            pass
        return out
    except Exception as e:
        print("Combine.Merge failed:", e)

    raise Exception("Union failed. Bodies may not intersect/touch, or are invalid.")


def cleanup_sheet_bodies(GetRootPart, Delete, BodySelection,
                         delete_all_non_solids=True, name_filter=None, verbose=True):
    root = GetRootPart()
    sheet_bodies = []

    for b in list(root.Bodies):
        try:
            is_solid = getattr(b, "IsSolid", None)
            is_sheet = getattr(b, "IsSheetBody", None)

            if is_sheet is None:
                is_sheet = (is_solid is False)

            if delete_all_non_solids:
                candidate = (is_sheet is True) or (is_solid is False)
            else:
                candidate = (is_solid is False)

            if candidate:
                if (name_filter is None) or (b.Name == name_filter):
                    sheet_bodies.append(b)
                    if verbose:
                        print("FOUND non-solid body:",
                              "Name='{}'".format(b.Name),
                              "IsSolid={}".format(is_solid),
                              "IsSheetBody={}".format(getattr(b, "IsSheetBody", "NA")))
        except Exception as e:
            if verbose:
                print("Skipping body due to exception:", e)

    if not sheet_bodies:
        if verbose:
            print("No sheet/surface bodies found to delete.")
        return

    Delete.Execute(BodySelection.Create(sheet_bodies))
    if verbose:
        print("Deleted {} sheet/surface bodies.".format(len(sheet_bodies)))

def get_or_create_component(GetRootPart, Window, Part, Component, name):
    root = GetRootPart()

    for c in list(root.Components):
        try:
            if c.GetName() == name:
                return c
        except:
            if getattr(c, "Name", "") == name:
                return c

    doc = Window.ActiveWindow.Document
    part_def = Part.Create(doc, name)
    comp = Component.Create(root, part_def)
    return comp


def move_bodies_to_component(ComponentHelper, BodySelection, bodies, target_component):
    bodies = [b for b in bodies if b is not None]
    if not bodies:
        return

    sel = BodySelection.Create(bodies)
    ComponentHelper.MoveBodiesToComponent(sel, target_component)


def find_body_by_name_anywhere(GetRootPart, name):
    root = GetRootPart()

    for b in list(root.Bodies):
        if getattr(b, "Name", "") == name:
            return b

    for comp in list(root.Components):
        try:
            for b in list(comp.Content.Bodies):
                if getattr(b, "Name", "") == name:
                    return b
        except:
            pass

    return None