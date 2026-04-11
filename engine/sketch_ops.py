import math


def MM(x):
    return x / 1000.0


def P2(Point2D, x, y):
    return Point2D.Create(MM(x), MM(y))


def clear_all_sketch_curves(GetRootPart, Delete, CurveSelection):
    root = GetRootPart()
    try:
        crvs = list(root.Curves)
        if crvs:
            Delete.Execute(CurveSelection.Create(crvs))
    except:
        pass


def set_sketch_plane_xy(ViewHelper, Plane):
    ViewHelper.SetSketchPlane(Plane.PlaneXY)


def set_sketch_plane_yz_at_x(Frame, Point, Direction, ViewHelper, Plane, x0):
    frame = Frame.Create(
        Point.Create(MM(x0), MM(0), MM(0)),
        Direction.DirY,
        Direction.DirZ
    )
    ViewHelper.SetSketchPlane(Plane.Create(frame))


def set_sketch_plane_xy_at_z(Frame, Point, Direction, ViewHelper, Plane, z0):
    frame = Frame.Create(
        Point.Create(MM(0), MM(0), MM(z0)),
        Direction.DirX,
        Direction.DirY
    )
    ViewHelper.SetSketchPlane(Plane.Create(frame))


def set_sketch_plane_zx_at_y(Frame, Point, Direction, ViewHelper, Plane, y0):
    frame = Frame.Create(
        Point.Create(MM(0), MM(y0), MM(0)),
        Direction.DirZ,
        Direction.DirX
    )
    ViewHelper.SetSketchPlane(Plane.Create(frame))


def solidify_sketch(ViewHelper, InteractionMode):
    ViewHelper.SetViewMode(InteractionMode.Solid)


def sketch_circle(SketchCircle, Point2D, cx, cy, r):
    return SketchCircle.Create(Point2D.Create(MM(cx), MM(cy)), MM(r)).CreatedCurves


def sketch_doubleD(SketchArc, SketchLine, Point2D, dx, R):
    c1 = SketchArc.CreateSweepArc(
        Point2D.Create(MM(dx), 0),
        Point2D.Create(MM(dx), MM(R)),
        Point2D.Create(MM(dx), MM(-R)),
        True
    ).CreatedCurves
    c2 = SketchArc.CreateSweepArc(
        Point2D.Create(MM(-dx), 0),
        Point2D.Create(MM(-dx), MM(R)),
        Point2D.Create(MM(-dx), MM(-R)),
        False
    ).CreatedCurves
    l1 = SketchLine.Create(
        Point2D.Create(MM(-dx), MM(R)),
        Point2D.Create(MM(dx), MM(R))
    ).CreatedCurves
    l2 = SketchLine.Create(
        Point2D.Create(MM(-dx), MM(-R)),
        Point2D.Create(MM(dx), MM(-R))
    ).CreatedCurves
    return list(c1) + list(c2) + list(l1) + list(l2)


def sketch_rectangle_xy(SketchLine, Point2D, xmin, xmax, ymin, ymax):
    SketchLine.Create(P2(Point2D, xmin, ymin), P2(Point2D, xmax, ymin))
    SketchLine.Create(P2(Point2D, xmax, ymin), P2(Point2D, xmax, ymax))
    SketchLine.Create(P2(Point2D, xmax, ymax), P2(Point2D, xmin, ymax))
    SketchLine.Create(P2(Point2D, xmin, ymax), P2(Point2D, xmin, ymin))


def sketch_rectangle_uv(SketchLine, Point2D, umin, umax, vmin, vmax):
    SketchLine.Create(P2(Point2D, umin, vmin), P2(Point2D, umax, vmin))
    SketchLine.Create(P2(Point2D, umax, vmin), P2(Point2D, umax, vmax))
    SketchLine.Create(P2(Point2D, umax, vmax), P2(Point2D, umin, vmax))
    SketchLine.Create(P2(Point2D, umin, vmax), P2(Point2D, umin, vmin))


def sketch_doubleD_shifted(SketchArc, SketchLine, Point2D, x0, y0, dx, R):
    SketchArc.CreateSweepArc(P2(Point2D, x0 + dx, y0), P2(Point2D, x0 + dx, y0 + R), P2(Point2D, x0 + dx, y0 - R), True)
    SketchArc.CreateSweepArc(P2(Point2D, x0 - dx, y0), P2(Point2D, x0 - dx, y0 + R), P2(Point2D, x0 - dx, y0 - R), False)
    SketchLine.Create(P2(Point2D, x0 - dx, y0 + R), P2(Point2D, x0 + dx, y0 + R))
    SketchLine.Create(P2(Point2D, x0 - dx, y0 - R), P2(Point2D, x0 + dx, y0 - R))


def sketch_doubleD_shifted_YZ(SketchArc, SketchLine, Point2D, y_center, z_center, dx, R):
    SketchArc.CreateSweepArc(P2(Point2D, y_center + dx, z_center),
                             P2(Point2D, y_center + dx, z_center + R),
                             P2(Point2D, y_center + dx, z_center - R),
                             True)
    SketchArc.CreateSweepArc(P2(Point2D, y_center - dx, z_center),
                             P2(Point2D, y_center - dx, z_center + R),
                             P2(Point2D, y_center - dx, z_center - R),
                             False)
    SketchLine.Create(P2(Point2D, y_center - dx, z_center + R), P2(Point2D, y_center + dx, z_center + R))
    SketchLine.Create(P2(Point2D, y_center - dx, z_center - R), P2(Point2D, y_center + dx, z_center - R))


def sketch_elliptic(SketchNurbs, SketchArc, Point2D, W_val, H_val, h, n_ellipse=80):
    D = W_val / 2.0
    H = H_val / 2.0

    C0 = D - H
    xc = C0 + h * H

    if not (0.0 < xc < D):
        raise Exception("Need 0 < xc < D. Got xc=%.6f, D=%.6f" % (xc, D))

    denom = (2.0 * D - xc)
    if abs(denom) < 1e-12:
        raise Exception("Invalid: 2D - xc ~= 0 causes division by zero in rs.")

    rs = (H**2 + (D - xc) * D) / denom

    if (xc - D + rs) <= 0:
        raise Exception("Invalid: xc - D + rs must be > 0 for rx.")

    rx = math.sqrt((xc * H**2) / (xc - D + rs))

    if rs <= 0 or rx <= 0:
        raise Exception("Invalid geometry: rs or rx <= 0 (rs=%.6f, rx=%.6f)" % (rs, rx))

    arg = 1.0 - (xc * xc) / (rx * rx)
    if arg <= 0:
        raise Exception("Ellipse does not reach x=xc (arg=%.6f). Check h/W/H." % arg)

    yc = H * math.sqrt(arg)
    xs = [(-xc + (2.0 * xc) * i / float(n_ellipse)) for i in range(n_ellipse + 1)]

    curves = []

    pts_top = [P2(Point2D, x, H * math.sqrt(max(0.0, 1.0 - (x*x)/(rx*rx)))) for x in xs]
    curves.extend(SketchNurbs.CreateFrom2DPoints(False, pts_top).CreatedCurves)

    pts_bot = [P2(Point2D, x, -H * math.sqrt(max(0.0, 1.0 - (x*x)/(rx*rx)))) for x in reversed(xs)]
    curves.extend(SketchNurbs.CreateFrom2DPoints(False, pts_bot).CreatedCurves)

    a = D - rs
    phi = math.atan2(yc, (xc - a))

    xr = a + rs * math.cos(phi)
    yr = rs * math.sin(phi)
    start_r = P2(Point2D, xr, +yr)
    end_r = P2(Point2D, xr, -yr)
    curves.extend(SketchArc.CreateSweepArc(P2(Point2D, +a, 0.0), start_r, end_r, True).CreatedCurves)

    xl = -a - rs * math.cos(phi)
    yl = rs * math.sin(phi)
    start_l = P2(Point2D, xl, +yl)
    end_l = P2(Point2D, xl, -yl)
    curves.extend(SketchArc.CreateSweepArc(P2(Point2D, -a, 0.0), start_l, end_l, False).CreatedCurves)

    return curves


def sketch_trapezoid(SketchLine, Point2D, cx, cy, angle_rad, h, w_inner, w_outer):
    ca = math.cos(angle_rad)
    sa = math.sin(angle_rad)

    ux, uy = ca, sa
    vx, vy = -sa, ca

    hi = 0.5 * h
    wi = 0.5 * w_inner
    wo = 0.5 * w_outer

    local = [
        (-hi, +wi),
        (+hi, +wo),
        (+hi, -wo),
        (-hi, -wi),
    ]

    pts = []
    for (u, v) in local:
        x = cx + u * ux + v * vx
        y = cy + u * uy + v * vy
        pts.append(P2(Point2D, x, y))

    SketchLine.Create(pts[0], pts[1])
    SketchLine.Create(pts[1], pts[2])
    SketchLine.Create(pts[2], pts[3])
    SketchLine.Create(pts[3], pts[0])


def sketch_profile(cfg, SketchNurbs, SketchArc, SketchLine, Point2D, dx, R, W_val, H_val):
    if cfg["is_elliptic"]:
        return sketch_elliptic(SketchNurbs, SketchArc, Point2D, W_val, H_val, cfg["h_mix"])
    return sketch_doubleD(SketchArc, SketchLine, Point2D, dx, R)


def sketch_wrap_loop_with_drains(cfg, SketchNurbs, SketchArc, SketchLine, Point2D, x_drain, offset):
    curves = []

    dx_shield = cfg["dx_shield"]
    R_shield = cfg["R_shield"]
    r_drain = cfg["r_drain"]
    h_mix = cfg["h_mix"]
    H_outer = cfg["H_outer"]
    W_outer = cfg["W_outer"]

    if not cfg["is_elliptic"]:
        d = abs(x_drain - dx_shield)
        R, r = R_shield, r_drain
        theta = math.asin((R - r) / d)
        Ro, ro = R + offset, r + offset
        tx_S, ty_S = Ro * math.sin(theta), Ro * math.cos(theta)
        tx_D, ty_D = ro * math.sin(theta), ro * math.cos(theta)

        curves.extend(SketchArc.CreateSweepArc(P2(Point2D, x_drain, 0), P2(Point2D, x_drain + tx_D, ty_D), P2(Point2D, x_drain + tx_D, -ty_D), True).CreatedCurves)
        curves.extend(SketchArc.CreateSweepArc(P2(Point2D, -x_drain, 0), P2(Point2D, -x_drain - tx_D, ty_D), P2(Point2D, -x_drain - tx_D, -ty_D), False).CreatedCurves)

        curves.extend(SketchLine.Create(P2(Point2D, x_drain + tx_D, ty_D), P2(Point2D, dx_shield + tx_S, ty_S)).CreatedCurves)
        curves.extend(SketchLine.Create(P2(Point2D, x_drain + tx_D, -ty_D), P2(Point2D, dx_shield + tx_S, -ty_S)).CreatedCurves)
        curves.extend(SketchLine.Create(P2(Point2D, -x_drain - tx_D, ty_D), P2(Point2D, -dx_shield - tx_S, ty_S)).CreatedCurves)
        curves.extend(SketchLine.Create(P2(Point2D, -x_drain - tx_D, -ty_D), P2(Point2D, -dx_shield - tx_S, -ty_S)).CreatedCurves)

        curves.extend(SketchArc.CreateSweepArc(P2(Point2D, dx_shield, 0), P2(Point2D, dx_shield + tx_S, ty_S), P2(Point2D, dx_shield, Ro), False).CreatedCurves)
        curves.extend(SketchArc.CreateSweepArc(P2(Point2D, dx_shield, 0), P2(Point2D, dx_shield + tx_S, -ty_S), P2(Point2D, dx_shield, -Ro), True).CreatedCurves)
        curves.extend(SketchArc.CreateSweepArc(P2(Point2D, -dx_shield, 0), P2(Point2D, -dx_shield, Ro), P2(Point2D, -dx_shield - tx_S, ty_S), False).CreatedCurves)
        curves.extend(SketchArc.CreateSweepArc(P2(Point2D, -dx_shield, 0), P2(Point2D, -dx_shield, -Ro), P2(Point2D, -dx_shield - tx_S, -ty_S), True).CreatedCurves)

        curves.extend(SketchLine.Create(P2(Point2D, -dx_shield, Ro), P2(Point2D, dx_shield, Ro)).CreatedCurves)
        curves.extend(SketchLine.Create(P2(Point2D, -dx_shield, -Ro), P2(Point2D, dx_shield, -Ro)).CreatedCurves)
        return curves

    ro = r_drain + offset
    Ho = (H_outer / 2.0) + offset
    Do = (W_outer / 2.0) + offset

    C0_o = Do - Ho
    xc_o = C0_o + h_mix * Ho
    rs_o = (Ho**2 + (Do - xc_o) * Do) / (2.0 * Do - xc_o)
    rx_o = math.sqrt((xc_o * Ho**2) / (xc_o - Do + rs_o))

    d = abs(x_drain)
    theta = math.asin((Ho - ro) / d)

    tx_D = ro * math.sin(theta)
    ty_D = ro * math.cos(theta)

    m = math.tan(theta)
    xt = (rx_o**2 * m) / math.sqrt(rx_o**2 * m**2 + Ho**2)
    yt = Ho * math.sqrt(max(0, 1.0 - (xt**2 / rx_o**2)))

    curves.extend(SketchArc.CreateSweepArc(P2(Point2D, x_drain, 0), P2(Point2D, x_drain + tx_D, ty_D), P2(Point2D, x_drain + tx_D, -ty_D), True).CreatedCurves)
    curves.extend(SketchArc.CreateSweepArc(P2(Point2D, -x_drain, 0), P2(Point2D, -x_drain - tx_D, ty_D), P2(Point2D, -x_drain - tx_D, -ty_D), False).CreatedCurves)

    curves.extend(SketchLine.Create(P2(Point2D, x_drain + tx_D, ty_D), P2(Point2D, xt, yt)).CreatedCurves)
    curves.extend(SketchLine.Create(P2(Point2D, x_drain + tx_D, -ty_D), P2(Point2D, xt, -yt)).CreatedCurves)
    curves.extend(SketchLine.Create(P2(Point2D, -x_drain - tx_D, ty_D), P2(Point2D, -xt, yt)).CreatedCurves)
    curves.extend(SketchLine.Create(P2(Point2D, -x_drain - tx_D, -ty_D), P2(Point2D, -xt, -yt)).CreatedCurves)

    n_step = 20
    xs = [(-xt + (2.0 * xt) * i / float(n_step)) for i in range(n_step + 1)]

    pts_top = [P2(Point2D, x, Ho * math.sqrt(max(0.0, 1.0 - (x*x)/(rx_o*rx_o)))) for x in xs]
    pts_bot = [P2(Point2D, x, -Ho * math.sqrt(max(0.0, 1.0 - (x*x)/(rx_o*rx_o)))) for x in reversed(xs)]

    curves.extend(SketchNurbs.CreateFrom2DPoints(False, pts_top).CreatedCurves)
    curves.extend(SketchNurbs.CreateFrom2DPoints(False, pts_bot).CreatedCurves)

    return curves