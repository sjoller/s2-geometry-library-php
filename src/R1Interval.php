<?php

/**
 *
 */
class R1Interval {
	/**
	 * @var
	 */
	private $lo;
	/**
	 * @var
	 */
	private $hi;

	/** Interval constructor. If lo > hi, the interval is empty. */
	public function __construct($lo, $hi) {
		$this->lo = $lo;
		$this->hi = $hi;
	}

	/**
	 * Returns an empty interval. (Any interval where lo > hi is considered
	 * empty.)
	 */
	public static function emptya(): R1Interval {
		return new R1Interval(1, 0);
	}

	/**
	 * Convenience method to construct an interval containing a single point.
	 */
	public static function fromPoint($p): R1Interval {
		return new R1Interval($p, $p);
	}

	/**
	 * Convenience method to construct the minimal interval containing the two
	 * given points. This is equivalent to starting with an empty interval and
	 * calling AddPoint() twice, but it is more efficient.
	 */
	public static function fromPointPair($p1, $p2): R1Interval {
		if ($p1 <= $p2) {
			return new R1Interval($p1, $p2);
		}

		return new R1Interval($p2, $p1);
	}

	/**
	 * @return mixed
	 */
	public function getLo() {
		return $this->lo;
	}

	/**
	 * @return mixed
	 */
	public function getHi() {
		return $this->hi;
	}

	/**
	 * Return true if the interval is empty, i.e. it contains no points.
	 */
	public function isEmpty(): bool {
		return $this->getLo() > $this->getHi();
	}

	/**
	 * Return the center of the interval. For empty intervals, the result is
	 * arbitrary.
	 */
	public function getCenter(): float {
		return 0.5 * ($this->getLo() + $this->getHi());
	}

	/**
	 * Return the length of the interval. The length of an empty interval is
	 * negative.
	 */
	public function getLength() {
		return $this->getHi() - $this->getLo();
	}

	/** Return true if this interval contains the interval 'y'.
	 *
	 * @param $p
	 * @return bool
	 */
	public function contains($p): bool {
		if ($p instanceof self) {
			$y = $p;
			if ($y->isEmpty()) {
				return true;
			}

			return $y->getLo() >= $this->getLo() && $y->getHi() <= $this->getHi();
		}

		return $p >= $this->getLo() && $p <= $this->getHi();
	}

	/**
	 * Return true if the interior of this interval contains the entire interval
	 * 'y' (including its boundary).
	 *
	 * @param $p
	 * @return bool
	 */
	public function interiorContains($p): bool {
		if ($p instanceof self) {
			$y = $p;
			if ($y->isEmpty()) {
				return true;
			}

			return $y->getLo() > $this->getLo() && $y->getHi() < $this->getHi();
		}

		return $p > $this->getLo() && $p < $this->getHi();
	}

	/**
	 * Return true if this interval intersects the given interval, i.e. if they
	 * have any points in common.
	 *
	 * @param R1Interval $y
	 * @return bool
	 */
	public function intersects(R1Interval $y): bool {
		if ($this->getLo() <= $y->getLo()) {
			return $y->getLo() <= $this->getHi() && $y->getLo() <= $y->getHi();
		}

		return $this->getLo() <= $y->getHi() && $this->getLo() <= $this->getHi();
	}

	/**
	 * Return true if the interior of this interval intersects any point of the
	 * given interval (including its boundary).
	 *
	 * @param R1Interval $y
	 * @return bool
	 */
	public function interiorIntersects(R1Interval $y): bool {
		return $y->getLo() < $this->getHi() && $this->getLo() < $y->getHi() && $this->getLo() < $this->getHi() && $y->getLo() <= $y->getHi();
	}

	/** Expand the interval so that it contains the given point "p".
	 *
	 * @param $p
	 * @return R1Interval
	 */
	public function addPoint($p): R1Interval {
		if ($this->isEmpty()) {
			return self::fromPoint($p);
		}

		if ($p < $this->getLo()) {
			return new R1Interval($p, $this->getHi());
		}

		if ($p > $this->getHi()) {
			return new R1Interval($this->getLo(), $p);
		}

		return new R1Interval($this->getLo(), $this->getHi());
	}

	/**
	 * Return an interval that contains all points with a distance "radius" of a
	 * point in this interval. Note that the expansion of an empty interval is
	 * always empty.
	 *
	 * @param $radius
	 * @return R1Interval
	 */
	public function expanded($radius): R1Interval {
		// assert (radius >= 0);
		if ($this->isEmpty()) {
			return $this;
		}

		return new R1Interval($this->getLo() - $radius, $this->getHi() + $radius);
	}

	/**
	 * Return the smallest interval that contains this interval and the given
	 * interval "y".
	 *
	 * @param R1Interval $y
	 * @return R1Interval
	 */
	public function union(R1Interval $y): R1Interval {
		if ($this->isEmpty()) {
			return $y;
		}

		if ($y->isEmpty()) {
			return $this;
		}

		return new R1Interval(min($this->getLo(), $y->getLo()), max($this->getHi(), $y->getHi()));
	}

	/**
	 * Return the intersection of this interval with the given interval. Empty
	 * intervals do not need to be special-cased.
	 *
	 * @param R1Interval $y
	 * @return R1Interval
	 */
	public function intersection(R1Interval $y): R1Interval {
		return new R1Interval(max($this->getLo(), $y->getLo()), min($this->getHi(), $y->getHi()));
	}

	/**
	 * @param $that
	 * @return bool
	 */
	public function equals($that): bool {
		if ($that instanceof self) {
			$y = $that;

			// Return true if two intervals contain the same set of points.
			return ($this->getLo() === $y->getLo() && $this->getHi() === $y->getHi()) || ($this->isEmpty() && $y->isEmpty());
		}

		return false;
	}

	/**
	 * @return int
	 */
	public function hashCode(): int {
		if ($this->isEmpty()) {
			return 17;
		}

		$value = 17;
		$value = 37 * $value + (float)$this->lo;
		$value = 37 * $value + (float)$this->hi;

		return $value ^ S2::unsignedRightShift($value, 32);
	}

	/**
	 * Return true if length of the symmetric difference between the two intervals
	 * is at most the given tolerance.
	 *
	 * @param R1Interval $y
	 * @param null       $maxError
	 * @return bool
	 */
	public function approxEquals(R1Interval $y, $maxError = null): bool {
		if ($maxError === null) {
			return $this->approxEquals($y, 1e-15);
		}
		if ($this->isEmpty()) {
			return $y->getLength() <= $maxError;
		}
		if ($y->isEmpty()) {
			return $this->getLength() <= $maxError;
		}

		return abs($y->getLo() - $this->getLo()) + abs($y->getHi() - $this->getHi()) <= $maxError;
	}

	/**
	 * @return string
	 */
	public function toString(): string {
		return "[" . $this->getLo() . ", " . $this->getHi() . "]";
	}
}
