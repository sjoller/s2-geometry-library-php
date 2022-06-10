<?php

/**
 * An abstract directed edge from one S2Point to another S2Point.
 */

class S2Edge {
    /** @var S2Point */
    private S2Point $start;
    /** @var S2Point */
    private S2Point $end;

    /**
     * @param S2Point $start
     * @param S2Point $end
     */
    public function __construct(S2Point $start, S2Point $end) {
        $this->start = $start;
        $this->end = $end;
    }

    /**
     * @return S2Point
     */
    public function getStart(): S2Point {
        return $this->start;
    }

	/**
	 * @return S2Point
	 */
	public function getEnd(): S2Point {
        return $this->end;
    }

	/**
	 * @return string
	 */
	public function toString(): string {
        return sprintf("Edge: (%s -> %s)\n   or [%s -> %s]", $this->start->toDegreesString(), $this->end->toDegreesString(), $this->start, $this->end );
    }

	/**
	 * @return void
	 */
	public function hashCode() {
        return $this->getStart()->hashCode() - $this->getEnd()->hashCode();
    }

	/**
	 * @param $o
	 * @return bool
	 */
	public function equals($o): bool {
        if (!($o instanceof self)) {
            return false;
        }
        $other = $o;
        return $this->getStart()->equals($other->getStart()) && $this->getEnd()->equals($other->getEnd());
    }
}
