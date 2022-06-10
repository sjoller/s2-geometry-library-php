<?php

class S2LatLngRect implements S2Region {
	/** @var R1Interval */
	private R1Interval $lat;
	/** @var S1Interval */
	private S1Interval $lng;

	/**
	 * Construct a rectangle from minimum and maximum latitudes and longitudes. If
	 * lo.lng() > hi.lng(), the rectangle spans the 180 degree longitude line.
	 */
	public function __construct($lo, $hi) {
		if ($lo instanceof S2LatLng && $hi instanceof S2LatLng) {
			$this->lat = new R1Interval($lo->lat()->radians(), $hi->lat()->radians());
			$this->lng = new S1Interval($lo->lng()->radians(), $hi->lng()->radians());
		}
		else if ($lo instanceof R1Interval && $hi instanceof S1Interval) {
			$this->lat = $lo;
			$this->lng = $hi;
		}
		// assert (isValid());
	}

	/** The canonical empty rectangle */
	public static function emptya(): S2LatLngRect {
		return new S2LatLngRect(R1Interval::emptya(), S1Interval::emptya());
	}

	/** The canonical full rectangle. */
	public static function full(): S2LatLngRect {
		return new S2LatLngRect(self::fullLat(), self::fullLng());
	}

	/** The full allowable range of latitudes. */
	public static function fullLat(): R1Interval {
		return new R1Interval(-S2::M_PI_2, S2::M_PI_2);
	}

	/**
	 * The full allowable range of longitudes.
	 */
	public static function fullLng(): S1Interval {
		return S1Interval::full();
	}

	/**
	 * Construct a rectangle from a center point (in lat-lng space) and size in
	 * each dimension. If size.lng() is greater than 360 degrees it is clamped,
	 * and latitudes greater than +/- 90 degrees are also clamped. So for example,
	 * FromCenterSize((80,170),(20,20)) -> (lo=(60,150),hi=(90,-170)).
	 *
	 * @param S2LatLng $center
	 * @param S2LatLng $size
	 * @return S2LatLngRect
	 */
	public static function fromCenterSize(S2LatLng $center, S2LatLng $size): S2LatLngRect {
		return self::fromPoint($center)->expanded($size->mul(0.5));
	}

	/** Convenience method to construct a rectangle containing a single point. */
	public static function fromPoint(S2LatLng $p): S2LatLngRect {
		// assert (p.isValid());
		return new S2LatLngRect($p, $p);
	}

	/**
	 * Convenience method to construct the minimal bounding rectangle containing
	 * the two given points. This is equivalent to starting with an empty
	 * rectangle and calling AddPoint() twice. Note that it is different than the
	 * S2LatLngRect(lo, hi) constructor, where the first point is always used as
	 * the lower-left corner of the resulting rectangle.
	 */
	public static function fromPointPair(S2LatLng $p1, S2LatLng $p2): S2LatLngRect {
		// assert (p1.isValid() && p2.isValid());
		return new S2LatLngRect(R1Interval::fromPointPair($p1->lat()->radians(), $p2->lat()->radians()), S1Interval::fromPointPair($p1->lng()->radians(), $p2->lng()->radians()));
	}

	/**
	 * Return a latitude-longitude rectangle that contains the edge from "a" to
	 * "b". Both points must be unit-length. Note that the bounding rectangle of
	 * an edge can be larger than the bounding rectangle of its endpoints.
	 */
	public static function fromEdge(S2Point $a, S2Point $b): S2LatLngRect {
		// assert (S2.isUnitLength(a) && S2.isUnitLength(b));
		$r = self::fromPointPair(new S2LatLng($a), new S2LatLng($b));

		// Check whether the min/max latitude occurs in the edge interior.
		// We find the normal to the plane containing AB, and then a vector "dir" in
		// this plane that also passes through the equator. We use RobustCrossProd
		// to ensure that the edge normal is accurate even when the two points are
		// very close together.
		$ab = S2::robustCrossProd($a, $b);
		$dir = S2Point::crossProd($ab, new S2Point(0, 0, 1));
		$da = $dir->dotProd($a);
		$db = $dir->dotProd($b);

		if ($da * $db >= 0) {
			// Minimum and maximum latitude are attained at the vertices.
			return $r;
		}

		// Minimum/maximum latitude occurs in the edge interior. This affects the
		// latitude bounds but not the longitude bounds.
		$absLat = acos(abs($ab->z / $ab->norm()));

		if ($da < 0) {
			return new S2LatLngRect(new R1Interval($r->lat()->getLo(), $absLat), $r->lng());
		}

		return new S2LatLngRect(new R1Interval(-$absLat, $r->lat()->getHi()), $r->lng());
	}

	/**
	 * Return true if the rectangle is valid, which essentially just means that
	 * the latitude bounds do not exceed Pi/2 in absolute value and the longitude
	 * bounds do not exceed Pi in absolute value.
	 *
	 */
	public function isValid(): bool {
		// The lat/lng ranges must either be both empty or both non-empty.
		return (abs($this->lat->getLo()) <= S2::M_PI_2 && abs($this->lat->getHi()) <= S2::M_PI_2 && $this->lng->isValid() && $this->lat->isEmpty() === $this->lng->isEmpty());
	}

	// Accessor methods.
	public function latLo(): S1Angle {
		return S1Angle::sradians($this->lat->getLo());
	}

	public function latHi(): S1Angle {
		return S1Angle::sradians($this->lat->getHi());
	}

	public function lngLo(): S1Angle {
		return S1Angle::sradians($this->lng->getLo());
	}

	public function lngHi(): S1Angle {
		return S1Angle::sradians($this->lng->getHi());
	}

	public function lat(): R1Interval {
		return $this->lat;
	}

	public function lng(): S1Interval {
		return $this->lng;
	}

	public function lo(): S2LatLng {
		return new S2LatLng($this->latLo(), $this->lngLo());
	}

	public function hi(): S2LatLng {
		return new S2LatLng($this->latHi(), $this->lngHi());
	}

	/**
	 * Return true if the rectangle is empty, i.e. it contains no points at all.
	 */
	public function isEmpty(): bool {
		return $this->lat->isEmpty();
	}

	/**
	 * Return true if the rectangle is full, i.e. it contains all points.
	 *
	 * @return bool
	 */
	public function isFull(): bool {
		return $this->lat->equals(self::fullLat()) && $this->lng->isFull();
	}

	/**
	 * Return true if lng_.lo() > lng_.hi(), i.e. the rectangle crosses the 180
	 * degree latitude line.
	 */
	public function isInverted(): bool {
		return $this->lng->isInverted();
	}

	/**
	 * Return the k-th vertex of the rectangle (k = 0,1,2,3) in CCW order.
	 *
	 * @param $k
	 * @throws RuntimeException
	 * @return S2LatLng
	 */
	public function getVertex($k): S2LatLng {
		//     Return the points in CCW order (SW, SE, NE, NW).
		switch ($k) {
			case 0:
				return S2LatLng::fromRadians($this->lat->getLo(), $this->lng->getLo());

			case 1:
				return S2LatLng::fromRadians($this->lat->getLo(), $this->lng->getHi());

			case 2:
				return S2LatLng::fromRadians($this->lat->getHi(), $this->lng->getHi());

			case 3:
				return S2LatLng::fromRadians($this->lat->getHi(), $this->lng->getLo());

			default:
				throw new RuntimeException("Invalid vertex index.");
		}
	}

	/**
	 * Return the center of the rectangle in latitude-longitude space (in general
	 * this is not the center of the region on the sphere).
	 */
	public function getCenter(): S2LatLng {
		return S2LatLng::fromRadians($this->lat->getCenter(), $this->lng->getCenter());
	}

	/**
	 * Return the minimum distance (measured along the surface of the sphere)
	 * from a given point to the rectangle (both its boundary and its interior).
	 * The latLng must be valid.
	 */
	public function getDistance(S2LatLng $p): S1Angle {
		// The algorithm here is the same as in getDistance(S2LagLngRect), only with simplified calculations.
		$a = $this;

		Preconditions . checkState(!$a->isEmpty());
		Preconditions . checkArgument($p->isValid());

		if ($a->lng()->contains($p->lng()->radians())) {
			return S1Angle::radians(max(0.0, max($p->lat()->radians() - $a->lat()->hi(), $a->lat()->lo() - $p->lat()->radians())));
		}

		$interval = new S1Interval($a->lng()->getHi(), $a->lng()->complement()->getCenter());
		$aLng = $a->lng()->getLo();
		if ($interval->contains($p->lng()->radians())) {
			$aLng = $a->lng()->getHi();
		}

		$lo = S2LatLng::fromRadians($a->lat()->getLo(), $aLng)->toPoint();
		$hi = S2LatLng::fromRadians($a->lat()->getHi(), $aLng)->toPoint();
		$loCrossHi = S2LatLng::fromRadians(0, $aLng - S2 . M_PI_2)->normalized()->toPoint();

		return S2EdgeUtil::getDistance($p->toPoint(), $lo, $hi, $loCrossHi);
	}

	/**
	 * Return the minimum distance (measured along the surface of the sphere) to
	 * the given S2LatLngRect. Both S2LatLngRects must be non-empty.
	 */
	public function getDistance(S2LatLngRect $other): S1Angle {
		$a = $this;
		$b = $other;

		Preconditions . checkState(!$a->isEmpty());
		Preconditions . checkArgument(!$b->isEmpty());

		// First, handle the trivial cases where the longitude intervals overlap.
		if ($a->lng()->intersects($b->lng())) {
			if ($a->lat()->intersects($b->lat())) {
				return S1Angle::radians(0);  // Intersection between a and b.
			}

			// We found an overlap in the longitude interval, but not in the latitude
			// interval. This means the shortest path travels along some line of
			// longitude connecting the high-latitude of the lower rect with the
			// low-latitude of the higher rect.
			if ($a->lat()->getLo() > $b->lat()->getHi()) {
				$lo = $b->latHi();
				$hi = $a->latLo();
			}
			else {
				$lo = $a->latHi();
				$hi = $b->latLo();
			}

			return S1Angle::radians($hi->radians() - $lo->radians());
		}

		// The longitude intervals don't overlap. In this case, the closest points
		// occur somewhere on the pair of longitudinal edges which are nearest in
		// longitude-space.
		$loHi = S1Interval::fromPointPair($a->lng()->getLo(), $b->lng()->getHi());
		$hiLo = S1Interval::fromPointPair($a->lng()->getHi(), $b->lng()->getLo());
		if ($loHi->getLength() < $hiLo->getLength()) {
			$aLng = $a->lngLo();
			$bLng = $b->lngHi();
		}
		else {
			$aLng = $a->lngHi();
			$bLng = $b->lngLo();
		}

// The shortest distance between the two longitudinal segments will include
// at least one segment endpoint. We could probably narrow this down further
// to a single point-edge distance by comparing the relative latitudes of the
// endpoints, but for the sake of clarity, we'll do all four point-edge
// distance tests.
		$aLo = new S2LatLng($a->latLo(), $aLng)->toPoint();
    $aHi = new S2LatLng($a->latHi(), $aLng)->toPoint();
    $aLoCrossHi = S2LatLng::fromRadians(0, $aLng->radians() - S2::M_PI_2)->normalized()->toPoint();
    $bLo = new S2LatLng($b->latLo(), $bLng)->toPoint();
    $bHi = new S2LatLng($b->latHi(), $bLng)->toPoint();
    $bLoCrossHi = S2LatLng::fromRadians(0, $bLng->radians() - S2::M_PI_2)->normalized()->toPoint();

    return S1Angle::min(S2EdgeUtil::getDistance($aLo, $bLo, $bHi, $bLoCrossHi), S1Angle::min(S2EdgeUtil::getDistance($aHi, $bLo, $bHi, $bLoCrossHi), S1Angle::min(S2EdgeUtil::getDistance($bLo, $aLo, $aHi, $aLoCrossHi), S2EdgeUtil::getDistance($bHi, $aLo, $aHi, $aLoCrossHi))));
  }

	/**
	 * Return the width and height of this rectangle in latitude-longitude space.
	 * Empty rectangles have a negative width and height.
	 */
	public function getSize(): S2LatLng {
		return S2LatLng::fromRadians($this->lat->getLength(), $this->lng->getLength());
	}

	/**
	 * More efficient version of Contains() that accepts a S2LatLng rather than an
	 * S2Point.
	 */
	public function contains($cell) {
		if ($cell instanceof S2LatLng) {
			return ($this->lat->contains($cell->lat()->radians()) && $this->lng->contains($cell->lng()->radians()));
		}

		if ($cell instanceof self) {
			return $this->lat->contains($cell->lat) && $this->lng->contains($cell->lng);
		}

		if ($cell instanceof S2Cell) {
			return $this->contains($cell->getRectBound());
		}

		if ($cell instanceof S2Point) {
			return self::contains(new S2LatLng($p));
		}
	}

	/**
	 * Return true if and only if the given point is contained in the interior of
	 * the region (i.e. the region excluding its boundary). The point 'p' does not
	 * need to be normalized.
	 */
//	public function interiorContains(S2Point $p): bool {
//		return self::interiorContains(new S2LatLng($p));
//	}

	/**
	 * More efficient version of InteriorContains() that accepts a S2LatLng rather
	 * than an S2Point.
	 */
//	public function interiorContains(S2LatLng $ll): bool {
//		// assert (ll.isValid());
//		return ($this->lat->interiorContains($ll->lat()->radians()) && $this->lng->interiorContains($ll->lng()->radians()));
//	}

	/**
	 * Return true if and only if the interior of this rectangle contains all
	 * points of the given other rectangle (including its boundary).
	 */
//	public function interiorContains(S2LatLngRect $other): bool {
//		return ($this->lat->interiorContains($other->lat) && $this->lng->interiorContains($other->lng));
//	}

	/**
	 * Return true if this rectangle and the given other rectangle have any points in common.
	 *
	 * @param S2LatLngRect $other
	 * @return bool
	 */
	public function intersects(S2LatLngRect $other): bool {
		$latInt = $this->lat->intersects($other->lat);
		$lngInt = $this->lng->intersects($other->lng);

		// echo var_export($latInt, true) . ' ' . var_export($lngInt, true) . "\n";
		return $latInt && $lngInt;
	}

	/**
	 * Returns true if this rectangle intersects the given cell. (This is an exact
	 * test and may be fairly expensive, see also MayIntersect below.)
	 */
	  public function intersects(S2Cell $cell): bool {
		// First we eliminate the cases where one region completely contains the
		// other. Once these are disposed of, then the regions will intersect
		// if and only if their boundaries intersect.

		if ($this->isEmpty()) {
		  return false;
		}
		if ($this->contains($cell->getCenter())) {
		  return true;
		}
		if ($cell->contains($this->getCenter()->toPoint())) {
		  return true;
		}

		// Quick rejection test (not required for correctness).
		if (!$this->intersects($cell->getRectBound())) {
		  return false;
		}

		// Now check whether the boundaries intersect. Unfortunately, a
		// latitude-longitude rectangle does not have straight edges -- two edges
		// are curved, and at least one of them is concave.

		// Precompute the cell vertices as points and latitude-longitudes.
		$cellV = new S2Point[4];
		$cellLl = new S2LatLng[4];
		for ($i = 0; $i < 4; ++$i) {
		  $cellV[$i] = $cell->getVertex($i); // Must be normalized.
		  $cellLl[$i] = new S2LatLng($cellV[$i]);
		  if ($this->contains($cellLl[$i])) {
			return true; // Quick acceptance test.
		  }
		}

		for ($i = 0; $i < 4; ++$i) {
		  $edgeLng = S1Interval::fromPointPair($cellLl[$i]->lng()->radians(), $cellLl[($i + 1) & 3]->lng()->radians());
		  if (!$this->lng->intersects($edgeLng)) {
			continue;
		  }

		  $a = $cellV[$i];
		  $b = $cellV[($i + 1) & 3];
		  if ($edgeLng->contains($this->lng->getLo())) {
			if ($this->intersectsLngEdge($a, $b, $this->lat, $this->lng->getLo())) {
			  return true;
			}
		  }
		  if ($edgeLng->contains($this->lng->getHi())) {
			if ($this->intersectsLngEdge($a, $b, $this->lat, $this->lng->getHi())) {
			  return true;
			}
		  }
		  if ($this->intersectsLatEdge($a, $b, $this->lat->getLo(), $this->lng)) {
			return true;
		  }
		  if ($this->intersectsLatEdge($a, $b, $this->lat->getHi(), $this->lng)) {
			return true;
		  }
		}
		return false;
	  }


	/**
	 * Return true if and only if the interior of this rectangle intersects any
	 * point (including the boundary) of the given other rectangle.
	 */
	public function interiorIntersects(S2LatLngRect $other): bool {
		return ($this->lat->interiorIntersects($other->lat) && $this->lng->interiorIntersects($other->lng));
	}

//	public function addPoint(S2Point $p): S2LatLngRect {
//		return self::addPoint(new S2LatLng($p));
//	}

	/**
	 * Increase the size of the bounding rectangle to include the given point. The rectangle is expanded by the minimum amount possible.
	 *
	 * @param S2LatLng $ll
	 * @return S2LatLngRect
	 */
	public function addPoint(S2LatLng $ll): S2LatLngRect {
		// assert (ll->isValid());
		$newLat = $this->lat->addPoint($ll->lat()->radians());
		$newLng = $this->lng->addPoint($ll->lng()->radians());

		return new S2LatLngRect($newLat, $newLng);
	}

	/**
	 * Return a rectangle that contains all points whose latitude distance from
	 * this rectangle is at most margin->lat(), and whose longitude distance from
	 * this rectangle is at most margin->lng(). In particular, latitudes are
	 * clamped while longitudes are wrapped. Note that any expansion of an empty
	 * interval remains empty, and both components of the given margin must be
	 * non-negative.
	 *
	 * NOTE: If you are trying to grow a rectangle by a certain *distance* on the
	 * sphere (e.g. 5km), use the ConvolveWithCap() method instead.
	 *
	 * @return S2LatLngRect
	 */
	public function expanded(S2LatLng $margin): S2LatLngRect {
		// assert (margin->lat()->radians() >= 0 && margin->lng()->radians() >= 0);
		if ($this->isEmpty()) {
			return $this;
		}

		return new S2LatLngRect(
			$this->lat->expanded($margin->lat()->radians())->intersection($this->fullLat()),
			$this->lng->expanded($margin->lng()->radians())
		);
	}

	/**
	 * Return the smallest rectangle containing the union of this rectangle and
	 * the given rectangle.
	 */
//  public S2LatLngRect union(S2LatLngRect other) {
//    return new S2LatLngRect(lat->union(other->lat), lng->union(other->lng));
//  }

	/**
	 * Return true if the latitude and longitude intervals of the two rectangles
	 * are the same up to the given tolerance (see r1interval->h and s1interval->h
	 * for details).
	 */
	public function approxEquals(S2LatLngRect $other, float $maxError): bool {
		return ($this->lat->approxEquals($other->lat, $maxError) && $this->lng->approxEquals($other->lng, $maxError));
	}

//	public function approxEquals(S2LatLngRect $other): bool {
//		return $this->approxEquals($other, 1e-15);
//	}

	//@Override
	public function hashCode(): int {
		$value = 17;
		$value = 37 * $value + $this->lat->hashCode();

		return (37 * $value + $this->lng->hashCode());
	}

	// //////////////////////////////////////////////////////////////////////
	// S2Region interface (see {@code S2Region} for details):

	//@Override
	public function clone(): S2Region {
		return new S2LatLngRect($this->lo(), $this->hi());
	}

	public function getCapBound(): S2Cap {
		// We consider two possible bounding caps, one whose axis passes
		// through the center of the lat-long rectangle and one whose axis
		// is the north or south pole. We return the smaller of the two caps.

		if ($this->isEmpty()) {
			echo __METHOD__ . " empty\n";

			return S2Cap::sempty();
		}

		$poleZ = null;
		$poleAngle = null;
		if ($this->lat->getLo() + $this->lat->getHi() < 0) {
			// South pole axis yields smaller cap.
			$poleZ = -1;
			$poleAngle = S2::M_PI_2 + $this->lat->getHi();
		}
		else {
			$poleZ = 1;
			$poleAngle = S2::M_PI_2 - $this->lat->getLo();
		}
		$poleCap = S2Cap::fromAxisAngle(new S2Point(0, 0, $poleZ), S1Angle::sradians($poleAngle));

		// For bounding rectangles that span 180 degrees or less in longitude, the
		// maximum cap size is achieved at one of the rectangle vertices. For
		// rectangles that are larger than 180 degrees, we punt and always return a
		// bounding cap centered at one of the two poles.
		$lngSpan = $this->lng->getHi() - $this->lng->getLo();
		if (S2::IEEEremainder($lngSpan, 2 * S2::M_PI) >= 0) {
			if ($lngSpan < 2 * S2::M_PI) {
				$midCap = S2Cap::fromAxisAngle($this->getCenter()->toPoint(), S1Angle::sradians(0));
				for ($k = 0; $k < 4; ++$k) {
					$midCap = $midCap->addPoint($this->getVertex($k)->toPoint());
				}
				if ($midCap->height() < $poleCap->height()) {
					return $midCap;
				}
			}
		}

		return $poleCap;
	}

	public
	function getRectBound(): S2LatLngRect {
		return $this;
	}

	/**
	 * This test is cheap but is NOT exact. Use Intersects() if you want a more
	 * accurate and more expensive test. Note that when this method is used by an
	 * S2RegionCoverer, the accuracy isn't all that important since if a cell may
	 * intersect the region then it is subdivided, and the accuracy of this method
	 * goes up as the cells get smaller.
	 */
	public
	function mayIntersect(S2Cell $cell): bool {
		// This test is cheap but is NOT exact (see s2latlngrect->h).

		$rb = $cell->getRectBound();

		// echo __METHOD__ . $cell . ' ' . $rb . "\n";

		return $this->intersects($rb);
	}

	/**
	 * Return true if the edge AB intersects the given edge of constant latitude.
	 */
	  private static boolean intersectsLatEdge(S2Point a, S2Point b, double lat,
	  S1Interval lng) {
	  // Return true if the segment AB intersects the given edge of constant
	  // latitude. Unfortunately, lines of constant latitude are curves on
	  // the sphere. They can intersect a straight edge in 0, 1, or 2 points.
	  // assert (S2->isUnitLength(a) && S2->isUnitLength(b));

	  // First, compute the normal to the plane AB that points vaguely north.
	  S2Point z = S2Point->normalize(S2->robustCrossProd(a, b));
	  if (z->z < 0) {
	  z = S2Point->neg(z);
	  }

	  // Extend this to an orthonormal frame (x,y,z) where x is the direction
	  // where the great circle through AB achieves its maximium latitude.
	  S2Point y = S2Point->normalize(S2->robustCrossProd(z, new S2Point(0, 0, 1)));
	  S2Point x = S2Point->crossProd(y, z);
	  // assert (S2->isUnitLength(x) && x->z >= 0);

	  // Compute the angle "theta" from the x-axis (in the x-y plane defined
	  // above) where the great circle intersects the given line of latitude.
	  double sinLat = Math->sin(lat);
	  if (Math->abs(sinLat) >= x->z) {
	  return false; // The great circle does not reach the given latitude.
	  }
	  // assert (x->z > 0);
	  double cosTheta = sinLat / x->z;
	  double sinTheta = Math->sqrt(1 - cosTheta * cosTheta);
	  double theta = Math->atan2(sinTheta, cosTheta);

	  // The candidate intersection points are located +/- theta in the x-y
	  // plane. For an intersection to be valid, we need to check that the
	  // intersection point is contained in the interior of the edge AB and
	  // also that it is contained within the given longitude interval "lng".

	  // Compute the range of theta values spanned by the edge AB.
	  S1Interval abTheta = S1Interval->fromPointPair(Math->atan2(
	  a->dotProd(y), a->dotProd(x)), Math->atan2(b->dotProd(y), b->dotProd(x)));

	  if (abTheta->contains(theta)) {
	  // Check if the intersection point is also in the given "lng" interval.
	  S2Point isect = S2Point->add(S2Point->mul(x, cosTheta), S2Point->mul(y,
	  sinTheta));
	  if (lng->contains(Math->atan2(isect->y, isect->x))) {
	  return true;
	  }
	  }
	  if (abTheta->contains(-theta)) {
	  // Check if the intersection point is also in the given "lng" interval.
	  S2Point intersection = S2Point->sub(S2Point->mul(x, cosTheta), S2Point->mul(y, sinTheta));
	  if (lng->contains(Math->atan2(intersection->y, intersection->x))) {
	  return true;
	  }
	  }
	  return false;

	  }

	public function __toString() {
		return sprintf("[Lo=%s, Hi=%s]", $this->lo(), $this->hi());
	}
}
