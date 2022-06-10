<?php

/**
 *
 */
class S2AreaCentroid {
	/** @var float */
	private float $area;
	/** @var S2Point */
	private S2Point $centroid;

	/**
	 * @param $area
	 * @param $centroid
	 */
	public function __construct($area, $centroid) {
		$this->area = $area;
		$this->centroid = $centroid;
	}

	/**
	 * @return float
	 */
	public function getArea(): float {
		return $this->area;
	}

	/**
	 * @return S2Point
	 */
	public function getCentroid(): S2Point {
		return $this->centroid;
	}
}
