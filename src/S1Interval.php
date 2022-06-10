<?php

class S1Interval {
	private $lo;
	private $hi;

	/**
	 * Both endpoints must be in the range -Pi to Pi inclusive. The value -Pi is
	 * converted internally to Pi except for the Full() and Empty() intervals.
	 *
	 * @param      $lo
	 * @param      $hi
	 * @param bool $checked
	 */
	public function __construct($lo, $hi = null, bool $checked = false) {
		if ($lo instanceof self) {
			$this->lo = $lo->lo;
			$this->hi = $lo->hi;
		}
		else {
			$newLo = $lo;
			$newHi = $hi;
			if ($checked === false) {
				if ($lo === -S2::M_PI && $hi !== S2::M_PI) {
					$newLo = S2::M_PI;
				}
				if ($hi === -S2::M_PI && $lo !== S2::M_PI) {
					$newHi = S2::M_PI;
				}
			}
			$this->lo = $newLo;
			$this->hi = $newHi;
		}
	}

	/**
	 * @return S1Interval
	 */
	public static function emptya(): S1Interval {
		return new S1Interval(S2::M_PI, -S2::M_PI, true);
	}

	/**
	 * @return S1Interval
	 */
	public static function full(): S1Interval {
		return new S1Interval(-S2::M_PI, S2::M_PI, true);
	}

	/**
	 * Convenience method to construct the minimal interval containing the two
	 * given points. This is equivalent to starting with an empty interval and
	 * calling AddPoint() twice, but it is more efficient.
	 *
	 * @param $p1
	 * @param $p2
	 * @return S1Interval
	 */
	public static function fromPointPair($p1, $p2): S1Interval {
		// assert (abs(p1) <= S2::M_PI && abs(p2) <= S2::M_PI);
		if ($p1 === -S2::M_PI) {
			$p1 = S2::M_PI;
		}

		if ($p2 === -S2::M_PI) {
			$p2 = S2::M_PI;
		}

		if (self::positiveDistance($p1, $p2) <= S2::M_PI) {
			return new S1Interval($p1, $p2, true);
		}

		return new S1Interval($p2, $p1, true);
	}

	public function getLo() {
		return $this->lo;
	}

	public function getHi() {
		return $this->hi;
	}

	/** Return true if the interval is empty, i.e. it contains no points. */
	public function isEmpty(): bool {
		return $this->getLo() - $this->getHi() === 2 * S2::M_PI;
	}

	/* Return true if lo() > hi(). (This is true for empty intervals.) */
	public function isInverted(): bool {
		return $this->getLo() > $this->getHi();
	}

	/**
	 * Return the midpoint of the interval. For full and empty intervals, the
	 * result is arbitrary.
	 */
	public function getCenter(): float {
		$center = 0.5 * ($this->getLo() + $this->getHi());

		if (!$this->isInverted()) {
			return $center;
		}

		// Return the center in the range (-Pi, Pi].
		return ($center <= 0) ? ($center + S2::M_PI) : ($center - S2::M_PI);
	}

	/**
	 * Return the length of the interval. The length of an empty interval is
	 * negative.
	 */
	public function getLength() {
		$length = $this->getHi() - $this->getLo();

		if ($length >= 0) {
			return $length;
		}

		$length += 2 * S2::M_PI;

		// Empty intervals have a negative length.
		return ($length > 0) ? $length : -1;
	}

	/** Return true if the interior of the interval contains the point 'p'. *#/
	 * public boolean interiorContains(double p) {
	 * // Works for empty, full, and singleton intervals.
	 * // assert (abs(p) <= S2::M_PI);
	 * if (p === -S2::M_PI) {
	 * p = S2::M_PI;
	 * }
	 *
	 * if (isInverted()) {
	 * return p > lo() || p < hi();
	 * } else {
	 * return (p > lo() && p < hi()) || isFull();
	 * }
	 * }
	 *
	 * /**
	 * Return true if the interval contains the given interval 'y'. Works for
	 * empty, full, and singleton intervals.
	 *#/
	 * public boolean contains(final S1Interval y) {
	 * // It might be helpful to compare the structure of these tests to
	 * // the simpler Contains(double) method above.
	 *
	 * if (isInverted()) {
	 * if (y.isInverted()) {
	 * return y.lo() >= lo() && y.hi() <= hi();
	 * }
	 * return (y.lo() >= lo() || y.hi() <= hi()) && !isEmpty();
	 * } else {
	 * if (y.isInverted()) {
	 * return isFull() || y.isEmpty();
	 * }
	 * return y.lo() >= lo() && y.hi() <= hi();
	 * }
	 * }
	 *
	 * /**
	 * Returns true if the interior of this interval contains the entire interval
	 * 'y'. Note that x.InteriorContains(x) is true only when x is the empty or
	 * full interval, and x.InteriorContains(S1Interval(p,p)) is equivalent to
	 * x.InteriorContains(p).
	 *#/
	 * public boolean interiorContains(final S1Interval y) {
	 * if (isInverted()) {
	 * if (!y.isInverted()) {
	 * return y.lo() > lo() || y.hi() < hi();
	 * }
	 * return (y.lo() > lo() && y.hi() < hi()) || y.isEmpty();
	 * } else {
	 * if (y.isInverted()) {
	 * return isFull() || y.isEmpty();
	 * }
	 * return (y.lo() > lo() && y.hi() < hi()) || isFull();
	 * }
	 * }
	 *
	 * /**
	 * Return true if the two intervals contain any points in common. Note that
	 * the point +/-Pi has two representations, so the intervals [-Pi,-3] and
	 * [2,Pi] intersect, for example.
	 */
	public function intersects(S1Interval $y): bool {
		if ($this->isEmpty() || $y->isEmpty()) {
			return false;
		}

		if ($this->isInverted()) {
			// Every non-empty inverted interval contains Pi.
			return $y->isInverted() || $y->getLo() <= $this->getHi() || $y->getHi() >= $this->getLo();
		}

		if ($y->isInverted()) {
			return $y->getLo() <= $this->getHi() || $y->getHi() >= $this->getLo();
		}

		return $y->getLo() <= $this->getHi() && $y->getHi() >= $this->getLo();
	}

	/**
	 * Return an interval that contains all points within a distance "radius" of
	 * a point in this interval. Note that the expansion of an empty interval is
	 * always empty. The radius must be non-negative.
	 */
	public function expanded($radius): S1Interval {
		// assert (radius >= 0);
		if ($this->isEmpty()) {
			return $this;
		}

		// Check whether this interval will be full after expansion, allowing
		// for a 1-bit rounding error when computing each endpoint.
		if ($this->getLength() + 2 * $radius >= 2 * S2::M_PI - 1e-15) {
			return self::full();
		}

		// NOTE(dbeaumont): Should this remainder be 2 * M_PI or just M_PI ??  double lo = IEEEremainder(lo() - radius, 2 * S2::M_PI);
		$lo = S2::IEEEremainder($this->getLo() - $radius, 2 * S2::M_PI);
		$hi = S2::IEEEremainder($this->getHi() + $radius, 2 * S2::M_PI);

		if ($lo === -S2::M_PI) {
			$lo = S2::M_PI;
		}

		return new S1Interval($lo, $hi);
	}

	/**
	 * Return true if two intervals contains the same set of points.
	 *
	 * @param Object $that
	 * @return bool
	 * @Override
	 */
	public function equals(Object $that): bool {
		if ($that instanceof self) {
			$thatInterval = $that;
			return $this->getLo() === $thatInterval->getLo() && $this->getHi() === $thatInterval->getHi();
		}

		return false;
	}

	/**
	 * @return int
	 * @Override
	 */
	public function hashCode(): int {
		$value = 17;
		$value = 37 * $value + (float)$this->getLo();
		$value = 37 * $value + (float)$this->getHi();
		return S2::unsignedRightShift($value,32) ^ $value;
	}

	/**
	 * @return string
	 * @Override
	 */
	public function toString(): string {
		return '[' . $this->getLo() . ', ' . $this->getHi() . ']';
	}

	/**
	 * Compute the distance from "a" to "b" in the range [0, 2*Pi). This is
	 * equivalent to (drem(b - a - S2::M_PI, 2 * S2::M_PI) + S2::M_PI), except that
	 * it is more numerically stable (it does not lose precision for very small
	 * positive distances).
	 *
	 * @param $a
	 * @param $b
	 * @return float|mixed
	 */
	public static function positiveDistance($a, $b) {
		$d = $b - $a;

		if ($d >= 0) {
			return $d;
		}

		// We want to ensure that if b === Pi and a === (-Pi + eps),
		// the return result is approximately 2*Pi and not zero.
		return ($b + S2::M_PI) - ($a - S2::M_PI);
	}
}
