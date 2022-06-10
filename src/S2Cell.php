<?php
/**
 *
 */
const MAX_CELL_SIZE = 1 << S2CellId::MAX_LEVEL;

/**
 *
 */
const MAX_ERROR = 1.0 / (1 << 51);

/**
 *
 */
define('POLE_MIN_LAT', asin(sqrt(1.0 / 3.0)) - MAX_ERROR);

/**
 *
 */
class S2Cell implements S2Region {
	/**
	 *
	 */
	public const MAX_CELL_SIZE = MAX_CELL_SIZE;

	/**
	 * @var
	 */
	private $face;

	/**
	 * @var
	 */
	private $level;

	/**
	 * @var
	 */
	private $orientation;

	/** @var S2CellId */
	private S2CellId $cellId;

	/**
	 * @var array
	 */
	private array $uv = [];

	/**
	 * Default constructor used only internally.
	 */
	public function __construct($p = null) {
		if ($p !== null) {
			switch(get_class($p)) {
				case 'S2Point':
					$this->init(S2CellId::fromPoint($p));
				break;
				case 'S2LatLng':
					$this->init(S2CellId::fromLatLng($p));
				break;
				case 'S2CellId':
					$this->init($p);
				break;
			}
		}
	}

	// This is a static method in order to provide named parameters.*/

	/**
	 * @param $face
	 * @param $pos
	 * @param $level
	 * @return S2Cell
	 */
	public static function fromFacePosLevel($face, $pos, $level): S2Cell {
		return new S2Cell(S2CellId::fromFacePosLevel($face, $pos, $level));
	}

	/**
	 * @return S2CellId
	 */
	public function getCellId(): S2CellId {
		return $this->cellId;
	}

	/**
	 * @return S2CellId
	 */
	public function id(): S2CellId {
		return $this->cellId;
	}

	/**
	 * @return mixed
	 */
	public function getFace() {
		return $this->face;
	}

	/**
	 * @return mixed
	 */
	public function getLevel() {
		return $this->level;
	}

	/**
	 * @return mixed
	 */
	public function getOrientation() {
		return $this->orientation;
	}

	/**
	 * @return bool
	 */
	public function isLeaf(): bool {
		return $this->level === S2CellId::MAX_LEVEL;
	}

	/**
	 * @param int $k
	 * @return S2Point
	 */
	public function getVertex(int $k): S2Point {
		return S2Point::normalize($this->getVertexRaw($k));
	}

	/**
	 * Return the k-th vertex of the cell (k = 0,1,2,3). Vertices are returned in
	 * CCW order. The points returned by GetVertexRaw are not necessarily unit
	 * length.
	 */
	public function getVertexRaw(int $k): S2Point {
		// Vertices are returned in the order SW, SE, NE, NW.
		return S2Projections::faceUvToXyz($this->face, $this->uv[0][($k >> 1) ^ ($k & 1)], $this->uv[1][$k >> 1]);
	}

	/**
	 * @param int $k
	 * @return S2Point
	 */
	public function getEdge(int $k): S2Point {
		return S2Point::normalize($this->getEdgeRaw($k));
	}

	/**
	 * @param int $k
	 * @return S2Point
	 */
	public function getEdgeRaw(int $k): S2Point {
		switch ($k) {
			case 0:
				return S2Projections::getVNorm($this->face, $this->uv[1][0]); // South
			case 1:
				return S2Projections::getUNorm($this->face, $this->uv[0][1]); // East
			case 2:
				return S2Point::neg(S2Projections::getVNorm($this->face, $this->uv[1][1])); // North
			default:
				return S2Point::neg(S2Projections::getUNorm($this->face, $this->uv[0][0])); // West
		}
	}

	/**
	 * Return the inward-facing normal of the great circle passing through the
	 * edge from vertex k to vertex k+1 (mod 4). The normals returned by
	 * GetEdgeRaw are not necessarily unit length.
	 *
	 *  If this is not a leaf cell, set children[0..3] to the four children of
	 * this cell (in traversal order) and return true. Otherwise returns false.
	 * This method is equivalent to the following:
	 *
	 *  for (pos=0, id=child_begin(); id !== child_end(); id = id.next(), ++pos)
	 * children[i] = S2Cell(id);
	 *
	 * except that it is more than two times faster.
	 *
	 * @param S2Cell[] $children
	 */
	public function subdivide(array &$children): bool {
		// This function is equivalent to just iterating over the child cell ids
		// and calling the S2Cell constructor, but it is about 2.5 times faster.

		if ($this->cellId->isLeaf()) {
			return false;
		}

		// Compute the cell midpoint in uv-space.
		$uvMid = $this->getCenterUV();

		// Create four children with the appropriate bounds.
		$id = $this->cellId->childBegin();
		for ($pos = 0; $pos < 4; ++$pos, $id = $id->next()) {
			$child = &$children[$pos];
			$child->face = $this->face;
			$child->level = $this->level + 1;
			$new_o = S2::posToOrientation($pos);
			$child->orientation = $this->orientation ^ $new_o;
			// echo "this-ori:" . $this->orientation . " new_o:" . $new_o . " res:" . $child->orientation . "\n";
			$child->cellId = $id;
			$ij = S2::posToIJ($this->orientation, $pos);
			for ($d = 0; $d < 2; ++$d) {
				// The dimension 0 index (i/u) is in bit 1 of ij.
				$m = 1 - (($ij >> (1 - $d)) & 1);
				$child->uv[$d][$m] = $uvMid->get($d);
				$child->uv[$d][1 - $m] = $this->uv[$d][1 - $m];
			}
		}

		return true;
	}

	/**
	 * Return the direction vector corresponding to the center in (s,t)-space of
	 * the given cell. This is the point at which the cell is divided into four
	 * subcells; it is not necessarily the centroid of the cell in (u,v)-space or
	 * (x,y,z)-space. The point returned by GetCenterRaw is not necessarily unit
	 * length.
	 *
	 * @return S2Point
	 */
	public function getCenter(): S2Point {
		return S2Point::normalize($this->getCenterRaw());
	}

	/**
	 * @return S2Point
	 */
	public function getCenterRaw(): S2Point {
		return $this->cellId->toPointRaw();
	}

	/**
	 * Return the center of the cell in (u,v) coordinates (see {@code
	 * S2Projections}). Note that the center of the cell is defined as the point
	 * at which it is recursively subdivided into four children; in general, it is
	 * not at the midpoint of the (u,v) rectangle covered by the cell
	 *
	 * @return R2Vector
	 */
	public function getCenterUV(): R2Vector {
		$i = 0;
		$j = 0;
		$null = null;
		$this->cellId->toFaceIJOrientation($i, $j, $null);
		$cellSize = 1 << (S2CellId::MAX_LEVEL - $this->level);

		// TODO(dbeaumont): Figure out a better naming of the variables here (and elsewhere).
		$si = ($i & -$cellSize) * 2 + $cellSize - self::MAX_CELL_SIZE;
		$x = S2Projections::stToUV((1.0 / self::MAX_CELL_SIZE) * $si);

		$sj = ($j & -$cellSize) * 2 + $cellSize - self::MAX_CELL_SIZE;
		$y = S2Projections::stToUV((1.0 / self::MAX_CELL_SIZE) * $sj);

		return new R2Vector($x, $y);
	}

	/**
	 * Return the average area for cells at the given level.
	 *
	 * @param int|null $level
	 * @return float
	 */
	public function averageArea(int $level = null): float {
		if ($level === null) {
			/**
			 * Return the average area of cells at this level. This is accurate to within
			 * a factor of 1.7 (for S2_QUADRATIC_PROJECTION) and is extremely cheap to
			 * compute.
			 */
			//$level = $this->getLevel();
			$level = $this->level;
		}

		return S2Projections::AVG_AREA->getValue($level);
	}

	/**
	 * Return the approximate area of this cell. This method is accurate to within
	 * 3% percent for all cell sizes and accurate to within 0.1% for cells at
	 * level 5 or higher (i.e. 300km square or smaller). It is moderately cheap to
	 * compute.
	 */
	public function approxArea(): float {
		// All cells at the first two levels have the same area.
		if ($this->level < 2) {
			return $this->averageArea($this->level);
		}

		// First, compute the approximate area of the cell when projected
		// perpendicular to its normal. The cross product of its diagonals gives
		// the normal, and the length of the normal is twice the projected area.
		$flatArea = 0.5 * S2Point::crossProd(S2Point::sub($this->getVertex(2), $this->getVertex(0)), S2Point::sub($this->getVertex(3), $this->getVertex(1)))->norm();

		// Now, compensate for the curvature of the cell surface by pretending
		// that the cell is shaped like a spherical cap. The ratio of the
		// area of a spherical cap to the area of its projected disc turns out
		// to be 2 / (1 + sqrt(1 - r*r)) where "r" is the radius of the disc.
		// For example, when r=0 the ratio is 1, and when r=1 the ratio is 2.
		// Here we set Pi*r*r === flat_area to find the equivalent disc.
		return $flatArea * 2 / (1 + sqrt(1 - min(S2::M_1_PI * $flatArea, 1.0)));
	}

	/**
	 * Return the area of this cell as accurately as possible. This method is more
	 * expensive, but it is accurate to 6 digits of precision even for leaf cells
	 * (whose area is approximately 1e-18).
	 *
	 * @return float
	 */
	public function exactArea(): float {
		$v0 = $this->getVertex(0);
		$v1 = $this->getVertex(1);
		$v2 = $this->getVertex(2);
		$v3 = $this->getVertex(3);

		return S2::area($v0, $v1, $v2) + S2::area($v0, $v2, $v3);
	}

	// //////////////////////////////////////////////////////////////////////
	// S2Region interface (see {@code S2Region} for details):

	/**
	 * @return S2Region
	 * @Override
	 */
	public function clone(): S2Region {
		$clone = new S2Cell();
		$clone->face = $this->face;
		$clone->level = $this->level;
		$clone->orientation = $this->orientation;
		$clone->uv = $this->uv->clone();

		return $clone;
	}

	/**
	 * @return S2Cap
	 */
	public function getCapBound(): S2Cap {
		// Use the cell center in (u,v)-space as the cap axis. This vector is
		// very close to GetCenter() and faster to compute. Neither one of these
		// vectors yields the bounding cap with minimal surface area, but they
		// are both pretty close.
		//
		// It's possible to show that the two vertices that are furthest from
		// the (u,v)-origin never determine the maximum cap size (this is a
		// possible future optimization).

		$u = 0.5 * ($this->uv[0][0] + $this->uv[0][1]);
		$v = 0.5 * ($this->uv[1][0] + $this->uv[1][1]);
		$cap = S2Cap::fromAxisHeight(S2Point::normalize(S2Projections::faceUvToXyz($this->face, $u, $v)), 0);

		for ($k = 0; $k < 4; ++$k) {
			$cap = $cap->addPoint($this->getVertex($k));
		}

		return $cap;
	}
	// We grow the bounds slightly to make sure that the bounding rectangle
	// also contains the normalized versions of the vertices. Note that the
	// maximum result magnitude is Pi, with a floating-point exponent of 1.
	// Therefore, adding or subtracting 2**-51 will always change the result.
	/**
	 *
	 */
	public const MAX_ERROR = MAX_ERROR;

	// The 4 cells around the equator extend to +/-45 degrees latitude at the
	// midpoints of their top and bottom edges. The two cells covering the
	// poles extend down to +/-35.26 degrees at their vertices.
	// adding kMaxError (as opposed to the C version) because of asin and atan2
	// round-off errors
	/**
	 *
	 */
	public const POLE_MIN_LAT = POLE_MIN_LAT;

	// 35.26 degrees

	/**
	 * @return S2LatLngRect
	 */
	public function getRectBound(): S2LatLngRect {
		if ($this->level > 0) {
			// Except for cells at level 0, the latitude and longitude extremes are
			// attained at the vertices. Furthermore, the latitude range is
			// determined by one pair of diagonally opposite vertices and the
			// longitude range is determined by the other pair.
			//
			// We first determine which corner (i,j) of the cell has the largest
			// absolute latitude. To maximize latitude, we want to find the point in
			// the cell that has the largest absolute z-coordinate and the smallest
			// absolute x- and y-coordinates. To do this we look at each coordinate
			// (u and v), and determine whether we want to minimize or maximize that
			// coordinate based on the axis direction and the cell's (u,v) quadrant.
			$u = $this->uv[0][0] + $this->uv[0][1];
			$v = $this->uv[1][0] + $this->uv[1][1];

			$i = $u > 0 ? 1 : 0;
			$j = $v > 0 ? 1 : 0;

			if (S2Projections::getUAxis($this->face)->z === 0) {
				$i = $u < 0 ? 1 : 0;
			}
			if (S2Projections::getVAxis($this->face)->z === 0) {
				$j = $v < 0 ? 1 : 0;
			}

			$lat = R1Interval::fromPointPair($this->getLatitude($i, $j), $this->getLatitude(1 - $i, 1 - $j));
			$lat = $lat->expanded(self::MAX_ERROR)->intersection(S2LatLngRect::fullLat());
			if ($lat->getLo() === -S2::M_PI_2 || $lat->getHi() === S2::M_PI_2) {
				return new S2LatLngRect($lat, S1Interval::full());
			}
			$lng = S1Interval::fromPointPair($this->getLongitude($i, 1 - $j), $this->getLongitude(1 - $i, $j));

			return new S2LatLngRect($lat, $lng->expanded(self::MAX_ERROR));
		}

		// The face centers are the +X, +Y, +Z, -X, -Y, -Z axes in that order.
		// assert (S2Projections.getNorm(face).get(face % 3) === ((face < 3) ? 1 : -1));
		switch ($this->face) {
			case 0:
				return new S2LatLngRect(
					new R1Interval(-S2::M_PI_4, S2::M_PI_4),
					new S1Interval(-S2::M_PI_4, S2::M_PI_4)
				);

			case 1:
				return new S2LatLngRect(
					new R1Interval(-S2::M_PI_4, S2::M_PI_4),
					new S1Interval(S2::M_PI_4, 3 * S2::M_PI_4)
				);

			case 2:
				return new S2LatLngRect(
					new R1Interval(POLE_MIN_LAT, S2::M_PI_2),
					new S1Interval(-S2::M_PI, S2::M_PI)
				);

			case 3:
				return new S2LatLngRect(
					new R1Interval(-S2::M_PI_4, S2::M_PI_4),
					new S1Interval(3 * S2::M_PI_4, -3 * S2::M_PI_4)
				);

			case 4:
				return new S2LatLngRect(
					new R1Interval(-S2::M_PI_4, S2::M_PI_4),
					new S1Interval(-3 * S2::M_PI_4, -S2::M_PI_4)
				);

			default:
				return new S2LatLngRect(
					new R1Interval(-S2::M_PI_2, -POLE_MIN_LAT),
					new S1Interval(-S2::M_PI, S2::M_PI)
				);
		}
	}

	/**
	 * @param S2Cell $cell
	 * @return mixed
	 */
	public function mayIntersect(S2Cell $cell) {
		return $this->cellId->intersects($cell->cellId);
	}

	/**
	 * @param $other
	 * @return bool|void
	 */
	public function contains($other) {
		// We can't just call XYZtoFaceUV, because for points that lie on the
		// boundary between two faces (i.e. u or v is +1/-1) we need to return
		// true for both adjacent cells.
		if ($other instanceof S2Point) {
			$uvPoint = S2Projections::faceXyzToUv($this->face, $other);
			if ($uvPoint === null) {
				return false;
			}

			return ($uvPoint->x() >= $this->uv[0][0] && $uvPoint->x() <= $this->uv[0][1] && $uvPoint->y() >= $this->uv[1][0] && $uvPoint->y() <= $this->uv[1][1]);
		}

		if ($other instanceof self) {
			return $this->cellId->contains($other->cellId);
		}
	}

	/**
	 * @param S2CellId $id
	 * @return void
	 */
	private function init(S2CellId $id): void {
		$this->cellId = $id;
		$ij = array(0, 0);
		$mOrientation = 0;

		// echo "   $mOrientation\n";
		$this->face = $id->toFaceIJOrientation($ij[0], $ij[1], $mOrientation);
		// echo ">> $mOrientation\n";
		$this->orientation = $mOrientation;
		$this->level = $id->level();
		$cellSize = 1 << (S2CellId::MAX_LEVEL - $this->level);
		for ($d = 0; $d < 2; ++$d) {
			// Compute the cell bounds in scaled (i,j) coordinates.
			$sijLo = ($ij[$d] & -$cellSize) * 2 - self::MAX_CELL_SIZE;
			$sijHi = $sijLo + $cellSize * 2;
			$this->uv[$d][0] = S2Projections::stToUV((1.0 / self::MAX_CELL_SIZE) * $sijLo);
			$this->uv[$d][1] = S2Projections::stToUV((1.0 / self::MAX_CELL_SIZE) * $sijHi);
		}
	}

	// Internal method that does the actual work in the constructors.

	/**
	 * @param $i
	 * @param $j
	 * @return float
	 */
	private function getLatitude($i, $j): float {
		$p = S2Projections::faceUvToXyz($this->face, $this->uv[0][$i], $this->uv[1][$j]);

		return atan2($p->z, sqrt($p->x * $p->x + $p->y * $p->y));
	}

	/**
	 * @param $i
	 * @param $j
	 * @return float
	 */
	private function getLongitude($i, $j): float {
		$p = S2Projections::faceUvToXyz($this->face, $this->uv[0][$i], $this->uv[1][$j]);

		return atan2($p->y, $p->x);
	}

	/**
	 * Return the latitude or longitude of the cell vertex given by (i,j),
	 * where 'i' and 'j' are either 0 or 1.
	 */
	public function __toString() {
		return sprintf('[%d, %d, %d, %s]', $this->face, $this->level, $this->orientation, $this->cellId);
	}

	/**
	 * @return int
	 * @Override
	 */
	public function hashCode(): int {
		$value = 17;
		$value = 37 * (37 * (37 * $value + $this->face) + $this->orientation) + $this->level;

		return 37 * $value + $this->getCellId()->hashCode();
	}

	/**
	 * @param Object $that
	 * @return bool
	 * @Override
	 */
	public function equals(object $that): bool {
		if ($that instanceof self) {
			return $this->face === $that->face && $this->level === $that->level && $this->orientation === $that->orientation && $this->cellId->equals($that->cellId);
		}

		return false;
	}

}
