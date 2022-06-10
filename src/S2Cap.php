<?php

const ROUND_UP = 1.0 + 1.0 / (1 << 52);

class S2Cap implements S2Region {
	/**
	 * Multiply a positive number by this constant to ensure that the result of a
	 * floating point operation is at least as large as the true
	 * infinite-precision result.
	 */
	public const ROUND_UP = ROUND_UP;

	private S2Point $axis;
	private $height;

	/**
	 * Caps may be constructed from either an axis and a height, or an axis and an angle. To avoid ambiguity, there are no public constructors
	 *
	 * @param S2Point|null $axis
	 * @param float|null   $height
	 */
	private function __construct(S2Point $axis = null, float $height = null) {
		if ($axis === null || $height === null) {
			$this->axis = new S2Point();
			$this->height = 0;
		}
		else {
			$this->axis = $axis;
			$this->height = $height;
		}
	}

	/**
	 * Create a cap given its axis and the cap height, i.e. the maximum projected
	 * distance along the cap axis from the cap center. 'axis' should be a
	 * unit-length vector.
	 *
	 * @param S2Point $axis
	 * @param float   $height
	 * @return S2Cap
	 */
	public static function fromAxisHeight(S2Point $axis, float $height): S2Cap {
		// assert (S2.isUnitLength(axis));
		return new S2Cap($axis, $height);
	}

	/**
	 * Create a cap given its axis and the cap opening angle, i.e. maximum angle
	 * between the axis and a point on the cap. 'axis' should be a unit-length
	 * vector, and 'angle' should be between 0 and 180 degrees.
	 *
	 * @param S2Point $axis
	 * @param S1Angle $angle
	 * @return S2Cap
	 */
	public static function fromAxisAngle(S2Point $axis, S1Angle $angle): S2Cap {
		// The height of the cap can be computed as 1-cos(angle), but this isn't
		// very accurate for angles close to zero (where cos(angle) is almost 1).
		// Computing it as 2*(sin(angle/2)**2) gives much better precision.

		// assert (S2.isUnitLength(axis));
		$d = sin(0.5 * $angle->radians());

		return new S2Cap($axis, 2 * $d * $d);
	}

	/**
	 * Create a cap given its axis and its area in steradians. 'axis' should be a
	 * unit-length vector, and 'area' should be between 0 and 4 * M_PI.
	 *
	 * @param S2Point $axis
	 * @param float   $area
	 * @return S2Cap
	 */
	public static function fromAxisArea(S2Point $axis, float $area): S2Cap {
		// assert (S2.isUnitLength(axis));
		return new S2Cap($axis, $area / (2 * S2::M_PI));
	}

	/** Return an empty cap, i.e. a cap that contains no points.
	 *
	 * @return S2Cap
	 */
	public static function empty(): S2Cap {
		return new S2Cap(new S2Point(1, 0, 0), -1);
	}

	/**
	 * Return a full cap, i.e. a cap that contains all points.
	 *
	 * @return S2Cap
	 */
	public static function full(): S2Cap {
		return new S2Cap(new S2Point(1, 0, 0), 2);
	}

	// Accessor methods.

	/**
	 * @return S2Point
	 */
	public function axis(): S2Point {
		return $this->axis;
	}

	/**
	 * @return float|int
	 */
	public function height() {
		return $this->height;
	}

	/**
	 * @return float|int
	 */
	public function area() {
		return 2 * S2::M_PI * max(0.0, $this->height);
	}

	/**
	 * Return the cap opening angle in radians, or a negative number for empty
	 * caps.
	 *
	 * @return S1Angle
	 */
	public function angle(): S1Angle {
		// This could also be computed as acos(1 - height_), but the following
		// formula is much more accurate when the cap height is small. It
		// follows from the relationship h = 1 - cos(theta) = 2 sin^2(theta/2).
		if ($this->isEmpty()) {
			return S1Angle::sradians(-1);
		}

		return S1Angle::sradians(2 * asin(sqrt(0.5 * $this->height)));
	}

	/**
	 * We allow negative heights (to represent empty caps) but not heights greater
	 * than 2.
	 */
	public function isValid(): bool {
		return S2::isUnitLength($this->axis) && $this->height <= 2;
	}

	/** Return true if the cap is empty, i.e. it contains no points.
	 *
	 * @return bool
	 */
	public function isEmpty(): bool {
		return $this->height < 0;
	}

	/**
	 * Return the complement of the interior of the cap. A cap and its complement
	 * have the same boundary but do not share any interior points. The complement
	 * operator is not a bijection, since the complement of a singleton cap
	 * (containing a single point) is the same as the complement of an empty cap.
	 *
	 * @return S2Cap
	 */
	public function complement(): S2Cap {
		// The complement of a full cap is an empty cap, not a singleton.
		// Also make sure that the complement of an empty cap has height 2.
		$cHeight = $this->isFull() ? -1 : 2 - max($this->height, 0.0);

		return self::fromAxisHeight(S2Point::neg($this->axis), $cHeight);
	}

	/**
	 * Return true if and only if this cap contains the given other cap (in a set
	 * containment sense, e.g. every cap contains the empty cap).
	 *
	 * @param S2Cap $other
	 * @return bool
	 */
	public function contains($other): bool {
		if ($this->isFull() || $other->isEmpty()) {
			return true;
		}

		if ($other instanceof S2LatLngRect) {
			// If the cap does not contain all cell vertices, return false.
			// We check the vertices before taking the Complement() because we can't
			// accurately represent the complement of a very small cap (a height
			// of 2-epsilon is rounded off to 2).
			$vertices = [];
			for ($k = 0; $k < 4; ++$k) {
				$vertices[$k] = $other->getVertex($k);
				if (!$this->contains($vertices[$k])) {
					return false;
				}
			}
			// Otherwise, return true if the complement of the cap does not intersect
			// the cell. (This test is slightly conservative, because technically we
			// want Complement().InteriorIntersects() here.)
			return !$this->complement()->intersects($other, $vertices);

		}

		return $this->angle()->radians() >= $this->axis->angle($other->axis) + $other->angle()->radians();
	}

	/**
	 * Return true if and only if the interior of this cap intersects the given
	 * other cap. (This relationship is not symmetric, since only the interior of
	 * this cap is used.)
	 *
	 * @param S2Cap $other
	 * @return bool
	 */
	public function interiorIntersects(S2Cap $other): bool {
		// Interior(X) intersects Y if and only if Complement(Interior(X))
		// does not contain Y.
		return !$this->complement()->contains($other);
	}

	/**
	 * Return true if and only if the given point is contained in the interior of
	 * the region (i.e. the region excluding its boundary). 'p' should be a
	 * unit-length vector.
	 *
	 * @param S2Point $p
	 * @return bool
	 */
	public function interiorContains(S2Point $p): bool {
		// assert (S2.isUnitLength(p));
		return $this->isFull() || S2Point::sub($this->axis, $p)->norm2() < 2 * $this->height;
	}

	/**
	 * Increase the cap height if necessary to include the given point. If the cap
	 * is empty the axis is set to the given point, but otherwise it is left
	 * unchanged. 'p' should be a unit-length vector.
	 *
	 * @param S2Point $p
	 * @return S2Cap
	 */
	public function addPoint(S2Point $p): S2Cap {
		// Compute the squared chord length, then convert it into a height.
		// assert (S2.isUnitLength(p));
		if ($this->isEmpty()) {
			return new S2Cap($p, 0);
		}

		// To make sure that the resulting cap actually includes this point,
		// we need to round up the distance calculation. That is, after
		// calling cap.AddPoint(p), cap.Contains(p) should be true.
		$dist2 = S2Point::sub($this->axis, $p)->norm2();
		$newHeight = max($this->height, self::ROUND_UP * 0.5 * $dist2);

		return new S2Cap($this->axis, $newHeight);
	}

	/**
	 * Increase the cap height if necessary to include "other". If the current cap is empty it is set to the given other cap.
	 *
	 * @param S2Cap $other
	 * @return S2Cap
	 */
	public function addCap(S2Cap $other): S2Cap {
		if ($this->isEmpty()) {
			return new S2Cap($other->axis, $other->height);
		}

		// See comments for FromAxisAngle() and AddPoint(). This could be
		// optimized by doing the calculation in terms of cap heights rather
		// than cap opening angles.
		$angle = $this->axis->angle($other->axis) + $other->angle()->radians();
		if ($angle >= S2::M_PI) {
			return new S2Cap($this->axis, 2); //Full cap
		}

		$d = sin(0.5 * $angle);
		$newHeight = max($this->height, ROUND_UP * 2 * $d * $d);

		return new S2Cap($this->axis, $newHeight);
	}

	// //////////////////////////////////////////////////////////////////////
	// S2Region interface (see {@code S2Region} for details):

	/**
	 * @return $this
	 */
	public function getCapBound(): S2Cap {
		return $this;
	}

	/**
	 * @return S2LatLngRect
	 */
	public function getRectBound(): S2LatLngRect {
		if ($this->isEmpty()) {
			return S2LatLngRect::emptya();
		}

		// Convert the axis to a (lat,lng) pair, and compute the cap angle.
		$axisLatLng = new S2LatLng($this->axis);
		$capAngle = $this->angle()->radians();

		$allLongitudes = false;
		$lat = array();
		$lng = array();
		$lng[0] = -S2::M_PI;
		$lng[1] = S2::M_PI;

		// Check whether cap includes the south pole.
		$lat[0] = $axisLatLng->lat()->radians() - $capAngle;
		if ($lat[0] <= -S2::M_PI_2) {
			$lat[0] = -S2::M_PI_2;
			$allLongitudes = true;
		}
		// Check whether cap includes the north pole.
		$lat[1] = $axisLatLng->lat()->radians() + $capAngle;
		if ($lat[1] >= S2::M_PI_2) {
			$lat[1] = S2::M_PI_2;
			$allLongitudes = true;
		}
		if (!$allLongitudes) {
			// Compute the range of longitudes covered by the cap. We use the law
			// of sines for spherical triangles. Consider the triangle ABC where
			// A is the north pole, B is the center of the cap, and C is the point
			// of tangency between the cap boundary and a line of longitude. Then
			// C is a right angle, and letting a,b,c denote the sides opposite A,B,C,
			// we have sin(a)/sin(A) = sin(c)/sin(C), or sin(A) = sin(a)/sin(c).
			// Here "a" is the cap angle, and "c" is the colatitude (90 degrees
			// minus the latitude). This formula also works for negative latitudes.
			//
			// The formula for sin(a) follows from the relationship h = 1 - cos(a).

			$sinA = sqrt($this->height * (2 - $this->height));
			$sinC = cos($axisLatLng->lat()->radians());
			if ($sinA <= $sinC) {
				$angleA = asin($sinA / $sinC);
				$lng[0] = S2::IEEEremainder($axisLatLng->lng()->radians() - $angleA, 2 * S2::M_PI);
				$lng[1] = S2::IEEEremainder($axisLatLng->lng()->radians() + $angleA, 2 * S2::M_PI);
			}
		}

		return new S2LatLngRect(
			new R1Interval($lat[0], $lat[1]),
			new S1Interval(
				$lng[0],
				$lng[1]
			)
		);
	}

	public function mayIntersect(S2Cell $cell): bool {
		// If the cap contains any cell vertex, return true.
		$vertices = array();
		for ($k = 0; $k < 4; ++$k) {
			$vertices[$k] = $cell->getVertex($k);
			if ($this->contains($vertices[$k])) {
				return true;
			}
		}

		return $this->intersects($cell, $vertices);
	}

	/**
	 * Return true if the cap axis and height differ by at most "max_error" from
	 * the given cap "other".
	 *
	 * @param S2Cap $other
	 * @param float $maxError
	 * @return bool
	 */
	public function approxEquals(S2Cap $other, float $maxError = 1e-14): bool {
		return ($this->axis->aequal($other->axis, $maxError) && abs($this->height - $other->height) <= $maxError) || ($this->isEmpty() && $other->height <= $maxError) || ($other->isEmpty() && $this->height <= $maxError) || ($this->isFull() && $other->height >= 2 - $maxError) || ($other->isFull() && $this->height >= 2 - $maxError);
	}

	/**
	 * @return string
	 * @Override
	 */
	public function toString(): string {
		  return '[Point = ' . $this->axis->toString() . ' Height = ' . $this->height . ']';
	  }
}
