import os
import json
import datetime


def load_input_file(path):
    if path is None:
        return {}

    if not os.path.exists(path):
        print("Input file not found, falling back to SpaceClaim Parameters: " + str(path))
        return {}

    with open(path, "r") as f:
        data = json.load(f)

    print("Loaded input file: " + path)
    return data


def make_param_getter(input_params, Parameters):
    def get_param(name, default_val):
        # 1) Prefer external flat input file
        if name in input_params:
            return input_params[name]

        # 2) If using a snapshot JSON, look inside "driving"
        try:
            if "driving" in input_params and name in input_params["driving"]:
                return input_params["driving"][name]
        except:
            pass

        # 3) Fall back to SpaceClaim Parameters object
        try:
            return float(Parameters[name])
        except:
            pass

        # 4) Fall back to attribute-style Parameters
        try:
            return float(getattr(Parameters, name))
        except:
            pass

        # 5) Default
        return default_val

    return get_param


def _resolve_bending_mode(opt):
    if opt == 0:
        return "good_3point"
    elif opt == 1:
        return "bad_3point"
    elif opt == 2:
        return "good_2pulleys"
    elif opt == 3:
        return "bad_2pulleys"
    raise Exception("Unknown bending_mode_option: %s" % opt)


def _safe_filename(base_name, ext):
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    return "{0}_{1}.{2}".format(base_name, timestamp, ext)


def _get_output_dir():
    out_dir = os.environ.get("TEMP", None)
    if not out_dir:
        out_dir = os.path.expanduser("~")
    return out_dir


def write_single_param_csv(param_dict, filename):
    try:
        path = os.path.join(_get_output_dir(), filename)
        with open(path, "w") as f:
            f.write("Name,Value\n")
            for k in sorted(param_dict.keys()):
                f.write("{0},{1}\n".format(k, param_dict[k]))
        print("CSV written to: " + path)
        return path
    except Exception as e:
        print("Failed to write CSV '{0}': {1}".format(filename, str(e)))
        return None


def write_param_json(driving, driven, metadata, filename):
    try:
        path = os.path.join(_get_output_dir(), filename)
        payload = {
            "driving": driving,
            "driven": driven,
            "metadata": metadata,
        }
        with open(path, "w") as f:
            json.dump(payload, f, indent=2, sort_keys=True)
        print("JSON written to: " + path)
        return path
    except Exception as e:
        print("Failed to write JSON '{0}': {1}".format(filename, str(e)))
        return None


def load_config(Parameters, input_file="input_case1.json", write_logs=True):
    input_params = load_input_file(input_file)
    get_param = make_param_getter(input_params, Parameters)

    cfg = {}
    cfg["input_file"] = input_file

    # top-level switches
    cfg["run_benchmark_rig"] = False
    cfg["run_stiffness"] = False

    # logic options
    cfg["filler_opt"] = int(get_param("filler_option", 0))
    cfg["drain_opt"] = int(get_param("has_drains", 0))
    cfg["is_a_doublet"] = int(get_param("is_a_doublet", 0))
    cfg["is_elliptic"] = int(get_param("is_elliptic", 1))
    print("DEBUG input_params keys =", sorted(input_params.keys()))
    print("DEBUG is_elliptic =", cfg["is_elliptic"])
    cfg["h_mix"] = float(get_param("mixing_factor", 0.7))
    cfg["bending_benchmark_opt"] = int(get_param("bending_mode_option", 0))

    # conductor / core
    cfg["D_cond"] = float(get_param("diam_cond", 1.0))
    cfg["core_overlap_input"] = float(get_param("core_overlap", 0.01))
    cfg["D_core"] = float(get_param("diam_core", 2.8))
    cfg["merge_cores"] = int(get_param("merge_cores", 0))

    # shield
    cfg["t_shield"] = float(get_param("t_shield", 0.1))
    cfg["W_outer"] = float(get_param("W_outer", 8.0))
    cfg["H_outer"] = float(get_param("H_outer", 5.0))

    # drains / overwrap
    cfg["D_drain"] = float(get_param("diam_drain", 1.0))
    cfg["t_overwrap"] = float(get_param("t_overwrap", 0.1))
    cfg["L_extrude"] = float(get_param("length_extrude", 50.0))

    # multilumen
    cfg["is_a_multilumen"] = int(get_param("is_a_multilumen", 0))
    cfg["multilumen_shape_opt"] = int(get_param("multilumen_shape", 0))
    cfg["multilumen_wall_thickness"] = float(get_param("multilumen_wall_thickness", 0.1))
    cfg["n_multilumen_cavity"] = int(get_param("n_multilumen_cavity", 9))

    # optional rig params
    cfg["Initial_Gap"] = None
    cfg["Nose_Diam"] = None
    cfg["Nose_Length"] = None

    # resolved logic
    cfg["bending_mode"] = _resolve_bending_mode(cfg["bending_benchmark_opt"])
    cfg["core_mode"] = "merged" if cfg["merge_cores"] else "separate"
    cfg["filler_mode"] = "shell" if cfg["filler_opt"] == 1 else "fill"
    cfg["multilumen_shape_mode"] = "trapezoidal" if cfg["multilumen_shape_opt"] == 1 else "circular"

    # shell policy: force touching cores
    core_overlap = cfg["core_overlap_input"]
    if cfg["filler_mode"] == "shell":
        core_overlap = 0.0

    if core_overlap < 0:
        raise Exception(
            "core_overlap must be >= 0. Negative values would create a real gap between cores."
        )

    min_overlap_fill = 0.001
    if cfg["filler_mode"] == "fill" and core_overlap < min_overlap_fill:
        raise Exception(
            "For filled second extrusion, cores need a positive overlap to avoid zero-thickness filler geometry. "
            "Got core_overlap=%.6f mm." % core_overlap
        )

    cfg["core_overlap_effective"] = core_overlap

    # derived values
    cfg["r_cond"] = cfg["D_cond"] / 2.0
    cfg["r_core"] = cfg["D_core"] / 2.0
    cfg["r_drain"] = cfg["D_drain"] / 2.0

    cfg["C2C"] = cfg["D_core"] - cfg["core_overlap_effective"]
    cfg["dx_cond"] = 0.5 * cfg["C2C"]

    cfg["W_in"] = cfg["W_outer"] - 2.0 * cfg["t_shield"]
    cfg["H_in"] = cfg["H_outer"] - 2.0 * cfg["t_shield"]
    cfg["R_in"] = cfg["H_in"] / 2.0
    cfg["dx_in"] = (cfg["W_in"] - cfg["H_in"]) / 2.0
    cfg["R_shield"] = cfg["H_outer"] / 2.0
    cfg["dx_shield"] = (cfg["W_outer"] - cfg["H_outer"]) / 2.0

    # optional rig defaults now that geometry exists
    cfg["Initial_Gap"] = float(get_param("Initial_Gap", 0.05 * (cfg["H_outer"] + 2.0 * cfg["t_overwrap"])))
    cfg["Nose_Diam"] = float(get_param("Nose_Diam", 1.5 * (cfg["H_outer"] + 2.0 * cfg["t_overwrap"])))
    cfg["Nose_Length"] = float(get_param("Nose_Length", 4.0 * (cfg["H_outer"] + 2.0 * cfg["t_overwrap"])))

    print("Driven C2C = %.6f mm" % cfg["C2C"])
    print("Driven dx_cond = %.6f mm" % cfg["dx_cond"])

    driving_params = {
        "diam_cond": cfg["D_cond"],
        "diam_core": cfg["D_core"],
        "core_overlap_input": cfg["core_overlap_input"],
        "diam_drain": cfg["D_drain"],
        "t_shield": cfg["t_shield"],
        "t_overwrap": cfg["t_overwrap"],
        "W_outer": cfg["W_outer"],
        "H_outer": cfg["H_outer"],
        "length_extrude": cfg["L_extrude"],
        "filler_option": cfg["filler_opt"],
        "has_drains": cfg["drain_opt"],
        "is_a_doublet": cfg["is_a_doublet"],
        "is_elliptic": cfg["is_elliptic"],
        "mixing_factor": cfg["h_mix"],
        "bending_mode_option": cfg["bending_benchmark_opt"],
        "merge_cores": cfg["merge_cores"],
        "is_a_multilumen": cfg["is_a_multilumen"],
        "multilumen_shape": cfg["multilumen_shape_opt"],
        "multilumen_wall_thickness": cfg["multilumen_wall_thickness"],
        "n_multilumen_cavity": cfg["n_multilumen_cavity"],
        "Initial_Gap": cfg["Initial_Gap"],
        "Nose_Diam": cfg["Nose_Diam"],
        "Nose_Length": cfg["Nose_Length"],
    }

    driven_params = {
        "filler_mode": cfg["filler_mode"],
        "core_mode": cfg["core_mode"],
        "bending_mode": cfg["bending_mode"],
        "multilumen_shape_mode": cfg["multilumen_shape_mode"],
        "core_overlap_effective": cfg["core_overlap_effective"],
        "r_cond": cfg["r_cond"],
        "r_core": cfg["r_core"],
        "r_drain": cfg["r_drain"],
        "C2C": cfg["C2C"],
        "dx_cond": cfg["dx_cond"],
        "W_in": cfg["W_in"],
        "H_in": cfg["H_in"],
        "R_in": cfg["R_in"],
        "dx_in": cfg["dx_in"],
        "R_shield": cfg["R_shield"],
        "dx_shield": cfg["dx_shield"],
    }

    metadata_params = {
        "timestamp": str(datetime.datetime.now()),
        "script_name": "spaceclaim_cable_builder",
        "input_file": input_file,
    }

    cfg["driving_params"] = driving_params
    cfg["driven_params"] = driven_params
    cfg["metadata_params"] = metadata_params

    if write_logs:
        write_single_param_csv(driving_params, _safe_filename("driving_params", "csv"))
        write_single_param_csv(driven_params, _safe_filename("driven_params", "csv"))
        write_param_json(driving_params, driven_params, metadata_params, _safe_filename("params_snapshot", "json"))

    return cfg