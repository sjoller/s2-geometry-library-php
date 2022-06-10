<?php

interface S2Region {
    /** Return a bounding spherical cap.
     * @return S2Cap
     */
    public function getCapBound(): S2Cap;

    /** Return a bounding latitude-longitude rectangle. */
    public function getRectBound();

	/**
	 * If this method returns true, the region completely contains the given cell.
	 * Otherwise, either the region does not contain the cell or the containment
	 * relationship could not be determined.
	 *
	 * @param S2Cell|S2Cap $other
	 */
    public function contains($other);

	/**
	 * If this method returns false, the region does not intersect the given cell.
	 * Otherwise, either region intersects the cell, or the intersection
	 * relationship could not be determined.
	 *
	 * @param S2Cell $cell
	 */
    public function mayIntersect(S2Cell $cell);
}
