from engine import sketch_ops, body_ops


def build_cable(api, cfg):
    """
    Build only the cable geometry for now.
    No named selections, no topology, no rigs, no stiffness.
    Returns a dict of bodies.
    """

    # -----------------------------
    # SpaceClaim API handles
    # -----------------------------
    GetRootPart = api["GetRootPart"]
    Delete = api["Delete"]
    BodySelection = api["BodySelection"]
    CurveSelection = api["CurveSelection"]
    Fill = api["Fill"]
    Selection = api["Selection"]
    ExtrudeFaceOptions = api["ExtrudeFaceOptions"]
    ExtrudeType = api["ExtrudeType"]
    ExtrudeFaces = api["ExtrudeFaces"]
    FaceSelection = api["FaceSelection"]
    Direction = api["Direction"]
    MM = api["MM"]
    ViewHelper = api["ViewHelper"]
    InteractionMode = api["InteractionMode"]
    Plane = api["Plane"]
    Frame = api["Frame"]
    Point = api["Point"]
    SketchCircle = api["SketchCircle"]
    SketchArc = api["SketchArc"]
    SketchLine = api["SketchLine"]
    SketchNurbs = api["SketchNurbs"]
    Point2D = api["Point2D"]
    Combine = api["Combine"]
    MakeSolidsOptions = api["MakeSolidsOptions"]

    bodies = {}

    # -----------------------------
    # Short aliases from cfg
    # -----------------------------
    D_cond = cfg["D_cond"]
    D_core = cfg["D_core"]
    D_drain = cfg["D_drain"]
    L_extrude = cfg["L_extrude"]

    r_cond = cfg["r_cond"]
    r_core = cfg["r_core"]
    r_drain = cfg["r_drain"]

    dx_cond = cfg["dx_cond"]

    W_outer = cfg["W_outer"]
    H_outer = cfg["H_outer"]
    t_shield = cfg["t_shield"]
    t_overwrap = cfg["t_overwrap"]

    W_in = cfg["W_in"]
    H_in = cfg["H_in"]
    R_in = cfg["R_in"]
    dx_in = cfg["dx_in"]

    R_shield = cfg["R_shield"]
    dx_shield = cfg["dx_shield"]

    filler_mode = cfg["filler_mode"]
    core_mode = cfg["core_mode"]
    is_a_doublet = cfg["is_a_doublet"]
    is_a_multilumen = cfg["is_a_multilumen"]
    multilumen_shape_opt = cfg["multilumen_shape_opt"]
    multilumen_wall_thickness = cfg["multilumen_wall_thickness"]
    n_multilumen_cavity = cfg["n_multilumen_cavity"]
    drain_opt = cfg["drain_opt"]

    # -----------------------------
    # local helpers
    # -----------------------------
    def clear_sketch():
        sketch_ops.clear_all_sketch_curves(GetRootPart, Delete, CurveSelection)

    def set_xy():
        sketch_ops.set_sketch_plane_xy(ViewHelper, Plane)

    def sketch_profile(dx, R, W_val, H_val):
        return sketch_ops.sketch_profile(
            cfg, SketchNurbs, SketchArc, SketchLine, Point2D,
            dx, R, W_val, H_val
        )

    def extrude_and_name(name, length_mm, pick_largest=True):
        return body_ops.extrude_and_name(
            GetRootPart, Delete, BodySelection, CurveSelection,
            Fill, Selection, ExtrudeFaceOptions, ExtrudeType, ExtrudeFaces,
            FaceSelection, Direction, MM,
            lambda: sketch_ops.solidify_sketch(ViewHelper, InteractionMode),
            clear_sketch,
            name, length_mm, pick_largest
        )

    def extrude_from_explicit_curves(name, length_mm, curves):
        return body_ops.extrude_from_explicit_curves(
            GetRootPart, Delete, BodySelection, CurveSelection,
            Fill, ExtrudeFaceOptions, ExtrudeType, ExtrudeFaces,
            FaceSelection, Direction, MM, ViewHelper,
            InteractionMode, clear_sketch,
            name, length_mm, curves
        )

    def find_body(name):
        return body_ops.find_body_by_name(GetRootPart, name)

    def union_bodies(name, body_list):
        return body_ops.union_bodies(Combine, Selection, name, body_list)

    # -----------------------------
    # (1) Conductors
    # -----------------------------
    c1 = find_body("conductor[1]")
    if c1 is None:
        set_xy()
        clear_sketch()
        sketch_ops.sketch_circle(SketchCircle, Point2D, -dx_cond, 0.0, r_cond)
        c1 = extrude_and_name("conductor[1]", L_extrude, True)
    bodies["conductor[1]"] = c1

    c2 = find_body("conductor[2]")
    if c2 is None:
        set_xy()
        clear_sketch()
        sketch_ops.sketch_circle(SketchCircle, Point2D, dx_cond, 0.0, r_cond)
        c2 = extrude_and_name("conductor[2]", L_extrude, True)
    bodies["conductor[2]"] = c2

    # -----------------------------
    # (2) Cores
    # -----------------------------
    def sketch_core_with_lumens(center_x):
        set_xy()
        clear_sketch()

        # outer core circle
        sketch_ops.sketch_circle(SketchCircle, Point2D, center_x, 0.0, r_core)

        # conductor hole
        sketch_ops.sketch_circle(SketchCircle, Point2D, center_x, 0.0, r_cond)

        if not is_a_multilumen:
            return

        if n_multilumen_cavity <= 0:
            return

        t_septum = multilumen_wall_thickness
        ri = r_cond + t_septum
        ro = r_core - t_septum

        if ro <= ri:
            return

        angle_step = (2.0 * 3.141592653589793) / float(n_multilumen_cavity)
        is_right = (center_x > 0.0)
        base_angle = 3.141592653589793 if is_right else 0.0

        hole_clear = 0.005  # mm safety margin

        for i in range(n_multilumen_cavity):
            phi = base_angle + (i * angle_step)

            if multilumen_shape_opt == 0:
                # circular lumens
                rp = 0.5 * (ri + ro)
                hx = center_x + rp * __import__("math").cos(phi)
                hy = rp * __import__("math").sin(phi)

                chord = 2.0 * rp * __import__("math").sin(angle_step / 2.0)
                d_hole = min(
                    ro - ri - hole_clear,
                    chord - t_septum - hole_clear
                )

                if d_hole > hole_clear:
                    sketch_ops.sketch_circle(SketchCircle, Point2D, hx, hy, 0.5 * d_hole)

            else:
                # trapezoidal / parallel-wall lumens
                d_half = 0.5 * t_septum
                if d_half >= ri or d_half >= ro:
                    continue

                alpha_i = __import__("math").asin(d_half / ri)
                alpha_o = __import__("math").asin(d_half / ro)

                p1_ang = phi + (angle_step / 2.0 - alpha_i)
                p4_ang = phi + (angle_step / 2.0 - alpha_o)
                p2_ang = phi - (angle_step / 2.0 - alpha_i)
                p3_ang = phi - (angle_step / 2.0 - alpha_o)

                p1 = sketch_ops.P2(Point2D, center_x + ri * __import__("math").cos(p1_ang), ri * __import__("math").sin(p1_ang))
                p2 = sketch_ops.P2(Point2D, center_x + ri * __import__("math").cos(p2_ang), ri * __import__("math").sin(p2_ang))
                p3 = sketch_ops.P2(Point2D, center_x + ro * __import__("math").cos(p3_ang), ro * __import__("math").sin(p3_ang))
                p4 = sketch_ops.P2(Point2D, center_x + ro * __import__("math").cos(p4_ang), ro * __import__("math").sin(p4_ang))

                SketchArc.CreateSweepArc(sketch_ops.P2(Point2D, center_x, 0.0), p1, p2, True)
                SketchLine.Create(p2, p3)
                SketchArc.CreateSweepArc(sketch_ops.P2(Point2D, center_x, 0.0), p3, p4, False)
                SketchLine.Create(p4, p1)

    if core_mode == "separate":
        score1 = find_body("single_core[1]")
        if score1 is None:
            sketch_core_with_lumens(-dx_cond)
            score1 = extrude_and_name("single_core[1]", L_extrude, True)
        bodies["single_core[1]"] = score1

        score2 = find_body("single_core[2]")
        if score2 is None:
            sketch_core_with_lumens(dx_cond)
            score2 = extrude_and_name("single_core[2]", L_extrude, True)
        bodies["single_core[2]"] = score2

    else:
        score = find_body("single_core_merged")
        if score is None:
            core1 = find_body("single_core[1]")
            core2 = find_body("single_core[2]")

            if core1 is None:
                sketch_core_with_lumens(-dx_cond)
                core1 = extrude_and_name("single_core[1]", L_extrude, True)

            if core2 is None:
                sketch_core_with_lumens(dx_cond)
                core2 = extrude_and_name("single_core[2]", L_extrude, True)

            score = union_bodies("single_core_merged", [core1, core2])

        bodies["single_core_merged"] = score

    # -----------------------------
    # (3) Second Extrusion
    # -----------------------------
    if is_a_doublet:
        second_extrusion = find_body("Second_Extrusion")
        if second_extrusion is None:
            set_xy()
            clear_sketch()
            sketch_profile(dx_in, R_in, W_in, H_in)
            sketch_ops.sketch_circle(SketchCircle, Point2D, -dx_cond, 0.0, r_cond)
            sketch_ops.sketch_circle(SketchCircle, Point2D, dx_cond, 0.0, r_cond)
            second_extrusion = extrude_and_name("Second_Extrusion", L_extrude, True)
        bodies["Second_Extrusion"] = second_extrusion

    elif filler_mode == "fill":
        second_extrusion = find_body("Second_Extrusion")
        if second_extrusion is None:
            set_xy()
            clear_sketch()
            sketch_profile(dx_in, R_in, W_in, H_in)
            sketch_ops.sketch_circle(SketchCircle, Point2D, -dx_cond, 0.0, r_core)
            sketch_ops.sketch_circle(SketchCircle, Point2D, dx_cond, 0.0, r_core)
            second_extrusion = extrude_and_name("Second_Extrusion", L_extrude, True)
        bodies["Second_Extrusion"] = second_extrusion

    elif filler_mode == "shell":
        second_extrusion = find_body("Second_Extrusion")
        second_extrusion_outer = find_body("Second_Extrusion_outer_tmp")
        second_extrusion_inner = find_body("Second_Extrusion_inner_tmp")

        # touching cores only for shell mode
        if abs(cfg["core_overlap_effective"]) > 1e-12:
            raise Exception("Shell mode expects touching cores only, but core_overlap_effective != 0.")

        # strict admissibility
        tol_shell_fit = 1e-6
        if (R_in - r_core) <= tol_shell_fit:
            raise Exception(
                "Shell mode impossible: no positive shell thickness available. "
                "Need R_in > r_core, but got R_in=%.6f mm and r_core=%.6f mm."
                % (R_in, r_core)
            )

        if abs(dx_cond - dx_in) > tol_shell_fit:
            suggested_W_outer = H_outer + D_core
            suggested_H_outer = W_outer - D_core
            raise Exception(
                "Shell mode geometry incompatible. Need dx_in = dx_cond for a uniform Double-D shell "
                "that also matches the shield envelope. "
                "Current dx_in=%.6f mm, dx_cond=%.6f mm. "
                "Equivalent requirement: W_outer - H_outer = D_core = %.6f mm. "
                "Suggested fixes: set W_outer=%.6f mm (keeping H_outer fixed) "
                "or set H_outer=%.6f mm (keeping W_outer fixed)."
                % (dx_in, dx_cond, D_core, suggested_W_outer, suggested_H_outer)
            )

        if second_extrusion is None:
            if second_extrusion_outer is None:
                set_xy()
                clear_sketch()
                outer_curves = sketch_ops.sketch_doubleD(SketchArc, SketchLine, Point2D, dx_cond, R_in)
                second_extrusion_outer = extrude_from_explicit_curves(
                    "Second_Extrusion_outer_tmp", L_extrude, outer_curves
                )

            if second_extrusion_inner is None:
                set_xy()
                clear_sketch()
                inner_curves = sketch_ops.sketch_doubleD(SketchArc, SketchLine, Point2D, dx_cond, r_core)
                second_extrusion_inner = extrude_from_explicit_curves(
                    "Second_Extrusion_inner_tmp", L_extrude, inner_curves
                )

            options = MakeSolidsOptions()
            options.SubtractFromTarget = True

            Combine.Intersect(
                BodySelection.Create([second_extrusion_outer]),
                BodySelection.Create([second_extrusion_inner]),
                options
            )

            try:
                Delete.Execute(BodySelection.Create([second_extrusion_inner]))
            except:
                pass

            second_extrusion = find_body("Second_Extrusion_outer_tmp")
            if second_extrusion is None:
                second_extrusion = second_extrusion_outer

            if second_extrusion is None:
                raise Exception("Could not refetch final Second_Extrusion body after boolean.")

            if second_extrusion is not None:
                second_extrusion.Name = "Second_Extrusion"

        bodies["Second_Extrusion"] = second_extrusion

    # -----------------------------
    # (4) Shield
    # -----------------------------
    shield = find_body("Shield")
    shield_outer = find_body("Shield_outer_tmp")
    shield_inner = find_body("Shield_inner_tmp")

    if shield is None:
        if shield_outer is None:
            set_xy()
            clear_sketch()
            outer_curves = sketch_profile(dx_shield, R_shield, W_outer, H_outer)
            shield_outer = extrude_from_explicit_curves("Shield_outer_tmp", L_extrude, outer_curves)

        if shield_inner is None:
            set_xy()
            clear_sketch()
            inner_curves = sketch_profile(dx_in, R_in, W_in, H_in)
            shield_inner = extrude_from_explicit_curves("Shield_inner_tmp", L_extrude, inner_curves)

        options = MakeSolidsOptions()
        options.SubtractFromTarget = True

        Combine.Intersect(
            BodySelection.Create([shield_outer]),
            BodySelection.Create([shield_inner]),
            options
        )

        try:
            Delete.Execute(BodySelection.Create([shield_inner]))
        except:
            pass

        shield = find_body("Shield_outer_tmp")
        if shield is None:
            shield = shield_outer

        if shield is not None:
            shield.Name = "Shield"
    bodies["Shield"] = shield

    # -----------------------------
    # (5) Drains
    # -----------------------------
    x_drain = dx_shield + R_shield + r_drain
    bodies["x_drain"] = x_drain

    if drain_opt:
        drain1 = find_body("drain[1]")
        if drain1 is None:
            set_xy()
            clear_sketch()
            sketch_ops.sketch_circle(SketchCircle, Point2D, -x_drain, 0.0, r_drain)
            drain1 = extrude_and_name("drain[1]", L_extrude, True)
        bodies["drain[1]"] = drain1

        drain2 = find_body("drain[2]")
        if drain2 is None:
            set_xy()
            clear_sketch()
            sketch_ops.sketch_circle(SketchCircle, Point2D, x_drain, 0.0, r_drain)
            drain2 = extrude_and_name("drain[2]", L_extrude, True)
        bodies["drain[2]"] = drain2

    # -----------------------------
    # (6) Overwrap
    # -----------------------------
    overwrap = find_body("Overwrap")
    overwrap_outer = find_body("Overwrap_outer_tmp")
    overwrap_inner = find_body("Overwrap_inner_tmp")

    if overwrap is None:
        if drain_opt:
            if overwrap_outer is None:
                set_xy()
                clear_sketch()
                outer_curves = sketch_ops.sketch_wrap_loop_with_drains(
                    cfg, SketchNurbs, SketchArc, SketchLine, Point2D, x_drain, t_overwrap
                )
                overwrap_outer = extrude_from_explicit_curves("Overwrap_outer_tmp", L_extrude, outer_curves)

            if overwrap_inner is None:
                set_xy()
                clear_sketch()
                inner_curves = sketch_ops.sketch_wrap_loop_with_drains(
                    cfg, SketchNurbs, SketchArc, SketchLine, Point2D, x_drain, 0.0
                )
                overwrap_inner = extrude_from_explicit_curves("Overwrap_inner_tmp", L_extrude, inner_curves)

        else:
            if overwrap_outer is None:
                set_xy()
                clear_sketch()
                outer_curves = sketch_profile(
                    dx_shield, R_shield + t_overwrap,
                    W_outer + 2.0 * t_overwrap,
                    H_outer + 2.0 * t_overwrap
                )
                overwrap_outer = extrude_from_explicit_curves("Overwrap_outer_tmp", L_extrude, outer_curves)

            if overwrap_inner is None:
                set_xy()
                clear_sketch()
                inner_curves = sketch_profile(dx_shield, R_shield, W_outer, H_outer)
                overwrap_inner = extrude_from_explicit_curves("Overwrap_inner_tmp", L_extrude, inner_curves)

        options = MakeSolidsOptions()
        options.SubtractFromTarget = True

        Combine.Intersect(
            BodySelection.Create([overwrap_outer]),
            BodySelection.Create([overwrap_inner]),
            options
        )

        try:
            Delete.Execute(BodySelection.Create([overwrap_inner]))
        except:
            pass

        overwrap = find_body("Overwrap_outer_tmp")
        if overwrap is None:
            overwrap = overwrap_outer

        if overwrap is not None:
            overwrap.Name = "Overwrap"


    bodies["Overwrap"] = overwrap

    return bodies