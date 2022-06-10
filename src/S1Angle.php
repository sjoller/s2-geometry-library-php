<?php

/**
 *
 */
class S1Angle {
	/** @var float */
	private $radians;

	/**
	 * @return float
	 */
	public function radians() {
		return $this->radians;
	}

	/**
	 * @param float $radians
	 * @return S1Angle
	 */
	public static function sradians(float $radians): S1Angle {
		return new S1Angle($radians);
	}

	/**
	 * @return float
	 */
	public function degrees() {
		return $this->radians * (180 / M_PI);
	}

	/**
	 * @param float $degrees
	 * @return S1Angle
	 */
	public static function sdegrees(float $degrees): S1Angle {
		return new S1Angle($degrees * (M_PI / 180));
	}

	/**
	 * @return float
	 */
	public function e5(): float {
		return round($this->degrees() * 1e5);
	}

	/**
	 * @return float
	 */
	public function e6(): float {
		return round($this->degrees() * 1e6);
	}

	/**
	 * @return float
	 */
	public function e7(): float {
		return round($this->degrees() * 1e7);
	}

	/**
	 * @param float|S2Point|null $radians_or_x
	 * @param S2Point|null        $y
	 * Return the angle between two points, which is also equal to the distance
	 * between these points on the unit sphere. The points do not need to be
	 * normalized.
	 */
	public function __construct($radians_or_x = null, S2Point $y = null) {
		if ($radians_or_x instanceof S2Point && $y instanceof S2Point) {
			$this->radians = $radians_or_x->angle($y);
		}
		else {
			//$this->radians = $radians_or_x === null ? 0 : $radians_or_x;
			$this->radians = $radians_or_x ?: 0;
		}
	}

	/**
	 * @param $that
	 * @return bool
	 */
	public function equals($that): bool {
		if ($that instanceof self) {
			return $this->radians() === $that->radians();
		}

		return false;
	}

	/**
	 * @return void
	 */
	public function hashCode(): void {
		//$value = Double.doubleToLongBits(radians);
		//return (int) (value ^ (value >>> 32));
	}

	/**
	 * @param S1Angle $that
	 * @return bool
	 */
	public function lessThan(S1Angle $that): bool {
		return $this->radians() < $that->radians();
	}

	/**
	 * @param S1Angle $that
	 * @return bool
	 */
	public function greaterThan(S1Angle $that): bool {
		return $this->radians() > $that->radians();
	}

	/**
	 * @param S1Angle $that
	 * @return bool
	 */
	public function lessOrEquals(S1Angle $that): bool {
		return $this->radians() <= $that->radians();
	}

	/**
	 * @param S1Angle $that
	 * @return bool
	 */
	public function greaterOrEquals(S1Angle $that): bool {
		return $this->radians() >= $that->radians();
	}

	/**
	 * @param S1Angle $left
	 * @param S1Angle $right
	 * @return S1Angle
	 */
	public static function max(S1Angle $left, S1Angle $right): S1Angle {
		return $right->greaterThan($left) ? $right : $left;
	}

	/**
	 * @param S1Angle $left
	 * @param S1Angle $right
	 * @return S1Angle
	 */
	public static function min(S1Angle $left, S1Angle $right): S1Angle {
		return $right->greaterThan($left) ? $left : $right;
	}

	/**
	 * @param $e5
	 * @return S1Angle
	 */
	public static function se5($e5): S1Angle {
		return self::sdegrees($e5 * 1e-5);
	}

	/**
	 * @param $e6
	 * @return S1Angle
	 */
	public static function se6($e6): S1Angle {
		// Multiplying by 1e-6 isn't quite as accurate as dividing by 1e6,
		// but it's about 10 times faster and more than accurate enough.
		return self::sdegrees($e6 * 1e-6);
	}

	/**
	 * @param $e7
	 * @return S1Angle
	 */
	public static function se7($e7): S1Angle {
		return self::sdegrees($e7 * 1e-7);
	}

	/**
	 * Writes the angle in degrees with a 'd' suffix, e.g. "17.3745d". By default
	 * 6 digits are printed; this can be changed using setprecision(). Up to 17
	 * digits are required to distinguish one angle from another.
	 */
	public function toString(): string {
		return $this->degrees() . 'd';
	}

	/**
	 * @param S1Angle $that
	 * @return int
	 */
	public function compareTo(S1Angle $that): int {
		//return $this->radians < $that->radians ? -1 : $this->radians > $that->radians ? 1 : 0;
		return $this->radians <=> $that->radians; // Spaceship operator available in php v7+
	}
}
