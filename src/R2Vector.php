<?php

/**
 *
 */
class R2Vector {
	/**
	 * @var int|mixed
	 */
	private $x;
	/**
	 * @var int|mixed
	 */
	private $y;

	/**
	 * @param $x
	 * @param $y
	 * @throws RuntimeException
	 */
	public function __construct($x = null, $y = null) {
		if ($x !== null && $y !== null) {
			$this->x = $x;
			$this->y = $y;
		}
		else if ($x !== null) {
			if (!is_array($x) || count($x) !== 2) {
				throw new RuntimeException("Points must have exactly 2 coordinates");
			}
			$this->x = $x[0];
			$this->y = $x[1];
		}
		else {
			$this->x = 0;
			$this->y = 0;
		}
	}

	/**
	 * @return int|mixed
	 */
	public function getX() {
		return $this->x;
	}

	/**
	 * @return int|mixed
	 */
	public function getY() {
		return $this->y;
	}

	/**
	 * @param $index
	 * @throws RuntimeException
	 * @return int|mixed
	 */
	public function get($index) {
		if ($index > 1) {
			throw new RuntimeException($index);
		}

		return $index === 0 ? $this->x : $this->y;
	}

	/**
	 * @param R2Vector $p1
	 * @param R2Vector $p2
	 * @throws RuntimeException
	 * @return R2Vector
	 */
	public static function add(R2Vector $p1, R2Vector $p2): R2Vector {
		return new R2Vector($p1->x + $p2->x, $p1->y + $p2->y);
	}

	/**
	 * @param R2Vector $p
	 * @param          $m
	 * @throws RuntimeException
	 * @return R2Vector
	 */
	public static function mul(R2Vector $p, $m): R2Vector {
		return new R2Vector($m * $p->x, $m * $p->y);
	}

	/**
	 * @return float|int
	 */
	public function norm2() {
		return ($this->x * $this->x) + ($this->y * $this->y);
	}

	/**
	 * @param R2Vector $p1
	 * @param R2Vector $p2
	 * @return float|int
	 */
	public static function sdotProd(R2Vector $p1, R2Vector $p2) {
		return ($p1->x * $p2->x) + ($p1->y * $p2->y);
	}

	/**
	 * @param R2Vector $that
	 * @return float|int
	 */
	public function dotProd(R2Vector $that) {
		return self::sdotProd($this, $that);
	}

	/**
	 * @param R2Vector $that
	 * @return float|int
	 */
	public function crossProd(R2Vector $that) {
		return $this->x * $that->y - $this->y * $that->x;
	}

	/**
	 * @param R2Vector $vb
	 * @return bool
	 */
	public function lessThan(R2Vector $vb): bool {
		if ($this->x < $vb->x) {
			return true;
		}
		if ($vb->x < $this->x) {
			return false;
		}
		if ($this->y < $vb->y) {
			return true;
		}

		return false;
	}

	/**
	 * @param $that
	 * @return bool
	 */
	public function equals($that): bool {
		if (!($that instanceof self)) {
			return false;
		}
		$thatPoint = $that;

		return $this->x === $thatPoint->x && $this->y === $thatPoint->y;
	}

	/**
	 * Calculates hashcode based on stored coordinates. Since we want +0.0 and
	 * -0.0 to be treated the same, we ignore the sign of the coordinates.
	 *
	 * @return int
	 */
	public function hashCode(): int {
		$value = 17;
		$value += 37 * $value + (float)abs($this->x);
		$value += 37 * $value + (float)abs($this->y);
		return $value ^ (S2::unsignedRightShift($value,32));
	}

	/**
	 * @return string
	 */
	public function toString(): string {
		return "(" . $this->x . ", " . $this->y . ")";
	}
}
