def _get_named_selection_by_name(GetRootPart, ns_name):
    root = GetRootPart()
    try:
        for ns in list(root.NamedSelections):
            try:
                if ns.GetName() == ns_name:
                    return ns
            except:
                if getattr(ns, "Name", "") == ns_name:
                    return ns
    except:
        pass
    return None


def _delete_named_selection_if_exists(GetRootPart, Delete, Selection, ns_name):
    ns = _get_named_selection_by_name(GetRootPart, ns_name)
    if ns is None:
        return
    try:
        Delete.Execute(Selection.Create(ns))
    except:
        try:
            ns.Delete()
        except:
            pass


def create_ns(GetRootPart, Delete, Selection, BodySelection, NamedSelection, name, body_or_list):
    if isinstance(body_or_list, list):
        items = body_or_list
    else:
        items = [body_or_list]

    items = [b for b in items if b is not None]
    if not items:
        print("create_ns skipped (no bodies) for:", name)
        return None

    ns_name = "NS_" + name

    _delete_named_selection_if_exists(GetRootPart, Delete, Selection, ns_name)

    sel = BodySelection.Create(items)
    res = NamedSelection.Create(sel, Selection.Empty())

    try:
        res.CreatedNamedSelection.SetName(ns_name)
    except:
        try:
            res.CreatedNamedSelection.Name = ns_name
        except Exception as e:
            print("WARNING: failed to rename Named Selection to", ns_name, ":", e)
            return None

    return res.CreatedNamedSelection


def create_all_body_named_selections(api, cfg, bodies):
    GetRootPart = api["GetRootPart"]
    Delete = api["Delete"]
    Selection = api["Selection"]
    BodySelection = api["BodySelection"]
    NamedSelection = api["NamedSelection"]

    created = {}

    created["NS_conductor[1]"] = create_ns(
        GetRootPart, Delete, Selection, BodySelection, NamedSelection,
        "conductor[1]", bodies.get("conductor[1]")
    )

    created["NS_conductor[2]"] = create_ns(
        GetRootPart, Delete, Selection, BodySelection, NamedSelection,
        "conductor[2]", bodies.get("conductor[2]")
    )

    if cfg["core_mode"] == "separate":
        created["NS_single_core[1]"] = create_ns(
            GetRootPart, Delete, Selection, BodySelection, NamedSelection,
            "single_core[1]", bodies.get("single_core[1]")
        )
        created["NS_single_core[2]"] = create_ns(
            GetRootPart, Delete, Selection, BodySelection, NamedSelection,
            "single_core[2]", bodies.get("single_core[2]")
        )
    else:
        created["NS_single_core_merged"] = create_ns(
            GetRootPart, Delete, Selection, BodySelection, NamedSelection,
            "single_core_merged", bodies.get("single_core_merged")
        )

    created["NS_Second_Extrusion"] = create_ns(
        GetRootPart, Delete, Selection, BodySelection, NamedSelection,
        "Second_Extrusion", bodies.get("Second_Extrusion")
    )

    created["NS_Shield"] = create_ns(
        GetRootPart, Delete, Selection, BodySelection, NamedSelection,
        "Shield", bodies.get("Shield")
    )

    if cfg["drain_opt"]:
        created["NS_drain[1]"] = create_ns(
            GetRootPart, Delete, Selection, BodySelection, NamedSelection,
            "drain[1]", bodies.get("drain[1]")
        )
        created["NS_drain[2]"] = create_ns(
            GetRootPart, Delete, Selection, BodySelection, NamedSelection,
            "drain[2]", bodies.get("drain[2]")
        )

    created["NS_Overwrap"] = create_ns(
        GetRootPart, Delete, Selection, BodySelection, NamedSelection,
        "Overwrap", bodies.get("Overwrap")
    )

    return created