<?php

class S2Projections {
    const S2_LINEAR_PROJECTION = 1;
    const S2_TAN_PROJECTION = 2;
    const S2_QUADRATIC_PROJECTION = 3;

    const S2_PROJECTION = self::S2_QUADRATIC_PROJECTION;

    // All of the values below were obtained by a combination of hand analysis and
    // Mathematica. In general, S2_TAN_PROJECTION produces the most uniform
    // shapes and sizes of cells, S2_LINEAR_PROJECTION is considerably worse, and
    // S2_QUADRATIC_PROJECTION is somewhere in between (but generally closer to
    // the tangent projection than the linear one).


    // The minimum area of any cell at level k is at least MIN_AREA.GetValue(k),
    // and the maximum is at most MAX_AREA.GetValue(k). The average area of all
    // cells at level k is exactly AVG_AREA.GetValue(k).
//  public static final Metric MIN_AREA = new Metric(2,
//    S2_PROJECTION === Projections.S2_LINEAR_PROJECTION ? 1 / (3 * Math.sqrt(3)) : // 0.192
//      S2_PROJECTION === Projections.S2_TAN_PROJECTION ? (S2.M_PI * S2.M_PI)
//        / (16 * S2.M_SQRT2) : // 0.436
//        S2_PROJECTION === Projections.S2_QUADRATIC_PROJECTION
//          ? 2 * S2.M_SQRT2 / 9 : // 0.314
//          0);
//  public static final Metric MAX_AREA = new Metric(2,
//    S2_PROJECTION === Projections.S2_LINEAR_PROJECTION ? 1 : // 1.000
//      S2_PROJECTION === Projections.S2_TAN_PROJECTION ? S2.M_PI * S2.M_PI / 16 : // 0.617
//        S2_PROJECTION === Projections.S2_QUADRATIC_PROJECTION
//          ? 0.65894981424079037 : // 0.659
//          0);
//  public static final Metric AVG_AREA = new Metric(2, S2.M_PI / 6); // 0.524)


    // Each cell is bounded by four planes passing through its four edges and
    // the center of the sphere. These metrics relate to the angle between each
    // pair of opposite bounding planes, or equivalently, between the planes
    // corresponding to two different s-values or two different t-values. For
    // example, the maximum angle between opposite bounding planes for a cell at
    // level k is MAX_ANGLE_SPAN.GetValue(k), and the average angle span for all
    // cells at level k is approximately AVG_ANGLE_SPAN.GetValue(k).
//  public static final Metric MIN_ANGLE_SPAN = new Metric(1,
//    S2_PROJECTION === Projections.S2_LINEAR_PROJECTION ? 0.5 : // 0.500
//      S2_PROJECTION === Projections.S2_TAN_PROJECTION ? S2.M_PI / 4 : // 0.785
//        S2_PROJECTION === Projections.S2_QUADRATIC_PROJECTION ? 2. / 3 : // 0.667
//          0);
//  public static final Metric MAX_ANGLE_SPAN = new Metric(1,
//    S2_PROJECTION === Projections.S2_LINEAR_PROJECTION ? 1 : // 1.000
//      S2_PROJECTION === Projections.S2_TAN_PROJECTION ? S2.M_PI / 4 : // 0.785
//        S2_PROJECTION === Projections.S2_QUADRATIC_PROJECTION
//          ? 0.85244858959960922 : // 0.852
//          0);
//  public static final Metric AVG_ANGLE_SPAN = new Metric(1, S2.M_PI / 4); // 0.785


    // The width of geometric figure is defined as the distance between two
    // parallel bounding lines in a given direction. For cells, the minimum
    // width is always attained between two opposite edges, and the maximum
    // width is attained between two opposite vertices. However, for our
    // purposes we redefine the width of a cell as the perpendicular distance
    // between a pair of opposite edges. A cell therefore has two widths, one
    // in each direction. The minimum width according to this definition agrees
    // with the classic geometric one, but the maximum width is different. (The
    // maximum geometric width corresponds to MAX_DIAG defined below.)
    //
    // For a cell at level k, the distance between opposite edges is at least
    // MIN_WIDTH.GetValue(k) and at most MAX_WIDTH.GetValue(k). The average
    // width in both directions for all cells at level k is approximately
    // AVG_WIDTH.GetValue(k).
    //
    // The width is useful for bounding the minimum or maximum distance from a
    // point on one edge of a cell to the closest point on the opposite edge.
    // For example, this is useful when "growing" regions by a fixed distance.
    public static function MIN_WIDTH() {
        if (self::S2_PROJECTION === self::S2_LINEAR_PROJECTION) $deriv = 1 / sqrt(6);
        else if (self::S2_PROJECTION === self::S2_TAN_PROJECTION) $deriv = S2::M_PI / (4 * S2::M_SQRT2);
        else if (self::S2_PROJECTION === self::S2_QUADRATIC_PROJECTION) $deriv = S2::M_SQRT2 / 3;
        else $deriv = 0;
        return new Metric(1, $deriv);
    }
//
//  public static final Metric MAX_WIDTH = new Metric(1, MAX_ANGLE_SPAN.deriv());
//  public static final Metric AVG_WIDTH = new Metric(1,
//    S2_PROJECTION === Projections.S2_LINEAR_PROJECTION ? 0.70572967292222848 : // 0.706
//      S2_PROJECTION === Projections.S2_TAN_PROJECTION ? 0.71865931946258044 : // 0.719
//        S2_PROJECTION === Projections.S2_QUADRATIC_PROJECTION
//          ? 0.71726183644304969 : // 0.717
//          0);

    // The minimum edge length of any cell at level k is at least
    // MIN_EDGE.GetValue(k), and the maximum is at most MAX_EDGE.GetValue(k).
    // The average edge length is approximately AVG_EDGE.GetValue(k).
    //
    // The edge length metrics can also be used to bound the minimum, maximum,
    // or average distance from the center of one cell to the center of one of
    // its edge neighbors. In particular, it can be used to bound the distance
    // between adjacent cell centers along the space-filling Hilbert curve for
    // cells at any given level.
//  public static final Metric MIN_EDGE = new Metric(1,
//    S2_PROJECTION === Projections.S2_LINEAR_PROJECTION ? S2.M_SQRT2 / 3 : // 0.471
//      S2_PROJECTION === Projections.S2_TAN_PROJECTION ? S2.M_PI / (4 * S2.M_SQRT2) : // 0.555
//        S2_PROJECTION === Projections.S2_QUADRATIC_PROJECTION ? S2.M_SQRT2 / 3 : // 0.471
//          0);
//  public static final Metric MAX_EDGE = new Metric(1, MAX_ANGLE_SPAN.deriv());
//  public static final Metric AVG_EDGE = new Metric(1,
//    S2_PROJECTION === Projections.S2_LINEAR_PROJECTION ? 0.72001709647780182 : // 0.720
//      S2_PROJECTION === Projections.S2_TAN_PROJECTION ? 0.73083351627336963 : // 0.731
//        S2_PROJECTION === Projections.S2_QUADRATIC_PROJECTION
//          ? 0.72960687319305303 : // 0.730
//          0);
    /*

  // The minimum diagonal length of any cell at level k is at least
  // MIN_DIAG.GetValue(k), and the maximum is at most MAX_DIAG.GetValue(k).
  // The average diagonal length is approximately AVG_DIAG.GetValue(k).
  //
  // The maximum diagonal also happens to be the maximum diameter of any cell,
  // and also the maximum geometric width (see the discussion above). So for
  // example, the distance from an arbitrary point to the closest cell center
  // at a given level is at most half the maximum diagonal length.
  public static final Metric MIN_DIAG = new Metric(1,
    S2_PROJECTION === Projections.S2_LINEAR_PROJECTION ? S2.M_SQRT2 / 3 : // 0.471
      S2_PROJECTION === Projections.S2_TAN_PROJECTION ? S2.M_PI / (3 * S2.M_SQRT2) : // 0.740
        S2_PROJECTION === Projections.S2_QUADRATIC_PROJECTION
          ? 4 * S2.M_SQRT2 / 9 : // 0.629
          0);
  public static final Metric MAX_DIAG = new Metric(1,
    S2_PROJECTION === Projections.S2_LINEAR_PROJECTION ? S2.M_SQRT2 : // 1.414
      S2_PROJECTION === Projections.S2_TAN_PROJECTION ? S2.M_PI / Math.sqrt(6) : // 1.283
        S2_PROJECTION === Projections.S2_QUADRATIC_PROJECTION
          ? 1.2193272972170106 : // 1.219
          0);
  public static final Metric AVG_DIAG = new Metric(1,
    S2_PROJECTION === Projections.S2_LINEAR_PROJECTION ? 1.0159089332094063 : // 1.016
      S2_PROJECTION === Projections.S2_TAN_PROJECTION ? 1.0318115985978178 : // 1.032
        S2_PROJECTION === Projections.S2_QUADRATIC_PROJECTION
          ? 1.03021136949923584 : // 1.030
          0);

  // This is the maximum edge aspect ratio over all cells at any level, where
  // the edge aspect ratio of a cell is defined as the ratio of its longest
  // edge length to its shortest edge length.
  public static final double MAX_EDGE_ASPECT =
      S2_PROJECTION === Projections.S2_LINEAR_PROJECTION ? S2.M_SQRT2 : // 1.414
      S2_PROJECTION === Projections.S2_TAN_PROJECTION ? S2.M_SQRT2 : // 1.414
      S2_PROJECTION === Projections.S2_QUADRATIC_PROJECTION ? 1.44261527445268292 : // 1.443
      0;

  // This is the maximum diagonal aspect ratio over all cells at any level,
  // where the diagonal aspect ratio of a cell is defined as the ratio of its
  // longest diagonal length to its shortest diagonal length.
  public static final double MAX_DIAG_ASPECT = Math.sqrt(3); // 1.732
*/
    public static function stToUV($s) {
        switch (self::S2_PROJECTION) {
            case self::S2_LINEAR_PROJECTION:
                return $s;

            case self::S2_TAN_PROJECTION:
                // Unfortunately, tan(M_PI_4) is slightly less than 1.0. This isn't due
                // to
                // a flaw in the implementation of tan(), it's because the derivative of
                // tan(x) at x=pi/4 is 2, and it happens that the two adjacent floating
                // point numbers on either side of the infinite-precision value of pi/4
                // have
                // tangents that are slightly below and slightly above 1.0 when rounded
                // to
                // the nearest double-precision result.
                $s = tan(S2::M_PI_4 * $s);
                return $s + (1.0 / (1 << 53)) * $s;

            case self::S2_QUADRATIC_PROJECTION:
                if ($s >= 0) {
                    return (1 / 3.) * ((1 + $s) * (1 + $s) - 1);
                } else {
                    return (1 / 3.) * (1 - (1 - $s) * (1 - $s));
                }
            default:
                throw new IllegalStateException("Invalid value for S2_PROJECTION");
        }
    }

    public static function uvToST($u) {
        switch (self::S2_PROJECTION) {
            case self::S2_LINEAR_PROJECTION:
                return $u;

            case self::S2_TAN_PROJECTION:
                return (4 * S2::M_1_PI) * atan($u);

            case self::S2_QUADRATIC_PROJECTION:
                if ($u >= 0) {
                    return sqrt(1 + 3 * $u) - 1;
                } else {
                    return 1 - sqrt(1 - 3 * $u);
                }
            default:
                throw new IllegalStateException("Invalid value for S2_PROJECTION");
        }
    }

    /**
     * Convert (face, u, v) coordinates to a direction vector (not necessarily
     * unit length).
     */
    public static function faceUvToXyz($face, $u, $v) {
        switch ($face) {
            case 0:
                return new S2Point(1, $u, $v);

            case 1:
                return new S2Point(-$u, 1, $v);

            case 2:
                return new S2Point(-$u, -$v, 1);

            case 3:
                return new S2Point(-1, -$v, -$u);

            case 4:
                return new S2Point($v, -1, -$u);

            default:
                return new S2Point($v, $u, -1);
        }
    }

    public static function validFaceXyzToUv($face, S2Point $p) {
        // assert (p.dotProd(faceUvToXyz(face, 0, 0)) > 0);
        switch ($face) {
            case 0:
                $pu = $p->y / $p->x;
                $pv = $p->z / $p->x;
                break;

            case 1:
                $pu = -$p->x / $p->y;
                $pv = $p->z / $p->y;
                break;

            case 2:
                $pu = -$p->x / $p->z;
                $pv = -$p->y / $p->z;
                break;

            case 3:
                $pu = $p->z / $p->x;
                $pv = $p->y / $p->x;
                break;

            case 4:
                $pu = $p->z / $p->y;
                $pv = -$p->x / $p->y;
                break;

            default:
                $pu = -$p->y / $p->z;
                $pv = -$p->x / $p->z;
                break;
        }
        return new R2Vector($pu, $pv);
    }

    public static function xyzToFace(S2Point $p) {
        $face = $p->largestAbsComponent();
        if ($p->get($face) < 0) {
            $face += 3;
        }
        return $face;
    }

    /*
  public static R2Vector faceXyzToUv(int face, S2Point p) {
    if (face < 3) {
      if (p.get(face) <= 0) {
        return null;
      }
    } else {
      if (p.get(face - 3) >= 0) {
        return null;
      }
    }
    return validFaceXyzToUv(face, p);
  }

  public static S2Point getUNorm(int face, double u) {
    switch (face) {
      case 0:
        return new S2Point(u, -1, 0);
      case 1:
        return new S2Point(1, u, 0);
      case 2:
        return new S2Point(1, 0, u);
      case 3:
        return new S2Point(-u, 0, 1);
      case 4:
        return new S2Point(0, -u, 1);
      default:
        return new S2Point(0, -1, -u);
    }
  }

  public static S2Point getVNorm(int face, double v) {
    switch (face) {
      case 0:
        return new S2Point(-v, 0, 1);
      case 1:
        return new S2Point(0, -v, 1);
      case 2:
        return new S2Point(0, -1, -v);
      case 3:
        return new S2Point(v, -1, 0);
      case 4:
        return new S2Point(1, v, 0);
      default:
        return new S2Point(1, 0, v);
    }
  }

  public static S2Point getNorm(int face) {
    return faceUvToXyz(face, 0, 0);
  }
*/
    public static function getUAxis($face) {
        switch ($face) {
            case 0:
                return new S2Point(0, 1, 0);

            case 1:
                return new S2Point(-1, 0, 0);

            case 2:
                return new S2Point(-1, 0, 0);

            case 3:
                return new S2Point(0, 0, -1);

            case 4:
                return new S2Point(0, 0, -1);

            default:
                return new S2Point(0, 1, 0);
        }
    }

    public static function getVAxis($face) {
        switch ($face) {
            case 0:
                return new S2Point(0, 0, 1);

            case 1:
                return new S2Point(0, 0, 1);

            case 2:
                return new S2Point(0, -1, 0);

            case 3:
                return new S2Point(0, -1, 0);

            case 4:
                return new S2Point(1, 0, 0);

            default:
                return new S2Point(1, 0, 0);
        }
    }
}
