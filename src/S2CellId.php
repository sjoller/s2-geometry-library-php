<?php

class S2CellId {

	// Although only 60 bits are needed to represent the index of a leaf
	// cell, we need an extra bit in order to represent the position of
	// the center of the leaf cell along the Hilbert curve.
	public const FACE_BITS = 3;
	public const NUM_FACES = 6;
	public const MAX_LEVEL = 30; // Valid levels: 0..MAX_LEVEL
	public const POS_BITS = 61; //2 * MAX_LEVEL + 1;
	public const MAX_SIZE = 0x40000000; //1 << MAX_LEVEL;

	// Constant related to unsigned long's
	public const MAX_UNSIGNED = -1; // Equivalent to 0xffffffffffffffffL

	// The following lookup tables are used to convert efficiently between an
	// (i,j) cell index and the corresponding position along the Hilbert curve.
	// "lookup_pos" maps 4 bits of "i", 4 bits of "j", and 2 bits representing the
	// orientation of the current cell into 8 bits representing the order in which
	// that subcell is visited by the Hilbert curve, plus 2 bits indicating the
	// new orientation of the Hilbert curve within that subcell. (Cell
	// orientations are represented as combination of kSwapMask and kInvertMask.)
	//
	// "lookup_ij" is an inverted table used for mapping in the opposite
	// direction.
	//
	// We also experimented with looking up 16 bits at a time (14 bits of position
	// plus 2 of orientation) but found that smaller lookup tables gave better
	// performance. (2KB fits easily in the primary cache.)

	// Values for these constants are *declared* in the *.h file. Even though
	// the declaration specifies a value for the constant, that declaration
	// is not a *definition* of storage for the value. Because the values are
	// supplied in the declaration, we don't need the values here. Failing to
	// define storage causes link errors for any code that tries to take the
	// address of one of these values.
	public const LOOKUP_BITS = 4;
	public const SWAP_MASK = 0x01;
	public const INVERT_MASK = 0x02;

	public static $LOOKUP_POS = null;
	public static $LOOKUP_IJ = null;

	/**
	 * This is the offset required to wrap around from the beginning of the
	 * Hilbert curve to the end or vice versa; see next_wrap() and prev_wrap().
	 */
	private static int $WRAP_OFFSET = (self::NUM_FACES) << self::POS_BITS;

	/**
	 * The id of the cell.
	 */
	public $id;

	public function __construct($id = null) {
		$this->id = $id ?? 0;
	}

	/** The default constructor returns an invalid cell id.
	 *
	 * @return S2CellId
	 */
	public static function none(): S2CellId {
		return new S2CellId();
	}

	/**
	 * Returns an invalid cell id guaranteed to be larger than any valid cell id.
	 * Useful for creating indexes.
	 *#/
	 * public static S2CellId sentinel() {
	 * return new S2CellId(MAX_UNSIGNED); // -1
	 * }
	 *
	 * /**
	 * Return a cell given its face (range 0..5), 61-bit Hilbert curve position
	 * within that face, and level (range 0..MAX_LEVEL). The given position will
	 * be modified to correspond to the Hilbert curve position at the center of
	 * the returned cell. This is a static function rather than a constructor in
	 * order to give names to the arguments.
	 *
	 * @param $face
	 * @param $pos
	 * @param $level
	 * @return S2CellId
	 */
	public static function fromFacePosLevel($face, $pos, $level): S2CellId {
		return (new S2CellId(($face << self::POS_BITS) + ($pos | 1)))->parent($level);
	}

	/**
	 * Return the leaf cell containing the given point (a direction vector, not
	 * necessarily unit length).
	 *
	 * @param S2Point $p
	 * @return S2CellId
	 */
	public static function fromPoint(S2Point $p): S2CellId {
		$face = S2Projections::xyzToFace($p);
		$uv = S2Projections::validFaceXyzToUv($face, $p);
		$i = self::stToIJ(S2Projections::uvToST($uv->getX()));
		$j = self::stToIJ(S2Projections::uvToST($uv->getY()));

		return self::fromFaceIJ($face, $i, $j);
	}

	/** Return the leaf cell containing the given S2LatLng. *#/
	 * public static S2CellId fromLatLng(S2LatLng ll) {
	 * return fromPoint(ll.toPoint());
	 * }
	 *
	 * public S2Point toPoint() {
	 * return S2Point.normalize(toPointRaw());
	 * }
	 *
	 * /**
	 * Return the direction vector corresponding to the center of the given cell.
	 * The vector returned by ToPointRaw is not necessarily unit length.
	 *
	 * @return S2Point
	 */
	public function toPointRaw(): S2Point {
		// First we compute the discrete (i,j) coordinates of a leaf cell contained
		// within the given cell. Given that cells are represented by the Hilbert
		// curve position corresponding at their center, it turns out that the cell
		// returned by ToFaceIJOrientation is always one of two leaf cells closest
		// to the center of the cell (unless the given cell is a leaf cell itself,
		// in which case there is only one possibility).
		//
		// Given a cell of size s >= 2 (i.e. not a leaf cell), and letting (imin,
		// jmin) be the coordinates of its lower left-hand corner, the leaf cell
		// returned by ToFaceIJOrientation() is either (imin + s/2, jmin + s/2)
		// (imin + s/2 - 1, jmin + s/2 - 1). We can distinguish these two cases by
		// looking at the low bit of "i" or "j". In the first case the low bit is
		// zero, unless s === 2 (i.e. the level just above leaf cells) in which case
		// the low bit is one.
		//
		// The following calculation converts (i,j) to the (si,ti) coordinates of
		// the cell center. (We need to multiply the coordinates by a factor of 2
		// so that the center of leaf cells can be represented exactly.)

		$i = 0;
		$j = 0;
		$null = null;
		$face = $this->toFaceIJOrientation($i, $j, $null);
		// System.out.println("i= " + i.intValue() + " j = " + j.intValue());
		if ($this->isLeaf()) {
			$delta = 1;
		}
		else {
			$delta = ((($i ^ ($this->id >> 2 & PHP_INT_MAX >> 1)) & 1) !== 0) ? 2 : 0;
			/* >>> */
		}
		$si = ($i << 1) + $delta - self::MAX_SIZE;
		$ti = ($j << 1) + $delta - self::MAX_SIZE;

		return self::faceSiTiToXYZ($face, $si, $ti);
	}

	/** Return the S2LatLng corresponding to the center of the given cell.
	 *
	 * @return S2LatLng
	 */
	public function toLatLng(): S2LatLng {
		return new S2LatLng($this->toPointRaw());
	}

	/**
	 * The 64-bit unique identifier for this cell.
	 */
	public function id() {
		return $this->id;
	}

	/**
	 * Return true if id() represents a valid cell.
	 *
	 * @return bool
	 */
	public function isValid(): bool {
		return $this->face() < self::NUM_FACES && (($this->lowestOnBit() & (float)'0x1555555555555555L') !== 0);
	}

	/**
	 * Which cube face this cell belongs to, in the range 0..5.
	 *
	 * @return int
	 */
	public function face(): int {
		return $this->id >> self::POS_BITS & PHP_INT_MAX >> (self::POS_BITS - 1);
		/* >>> */
	}

	/**
	 * The position of the cell center along the Hilbert curve over this face, in
	 * the range 0..(2**kPosBits-1).
	 *
	 * @return int
	 */
	public function pos(): int {
		return $this->id & (-1 >> self::FACE_BITS) & (PHP_INT_MAX >> (self::FACE_BITS - 1));
		/* >>> logical shift right */
	}

	/**
	 * Return the subdivision level of the cell (range 0..MAX_LEVEL).
	 *
	 * @return int
	 */
	public function level(): int {
		// Fast path for leaf cells.
		if ($this->isLeaf()) {
			return self::MAX_LEVEL;
		}
		$x = $this->id & 0xffffffff;
		$level = -1;
		if ($x !== 0) {
			$level += 16;
		}
		else {
			$x = $this->id >> 32 & PHP_INT_MAX >> 31;
			/* >>> */
		}
		// We only need to look at even-numbered bits to determine the
		// level of a valid cell id.
		$x &= -$x; // Get the lowest bit.
		if (($x & 0x00005555) !== 0) {
			$level += 8;
		}
		if (($x & 0x00550055) !== 0) {
			$level += 4;
		}
		if (($x & 0x05050505) !== 0) {
			$level += 2;
		}
		if (($x & 0x11111111) !== 0) {
			++$level;
		}

		// assert (level >= 0 && level <= MAX_LEVEL);
		return $level;
	}

	/**
	 * Return true if this is a leaf cell (more efficient than checking whether
	 * level() === MAX_LEVEL).
	 *
	 * @return bool
	 */
	public function isLeaf(): bool {
		return ($this->id & 1) !== 0;
	}

	/**
	 * Return true if this is a top-level face cell (more efficient than checking
	 * whether level() === 0).
	 *
	 * @return bool
	 */
	public function isFace(): bool {
		return ($this->id & (self::lowestOnBitForLevel(0) - 1)) === 0;
	}

	/**
	 * Return the child position (0..3) of this cell's ancestor at the given
	 * level, relative to its parent. The argument should be in the range
	 * 1..MAX_LEVEL. For example, child_position(1) returns the position of this
	 * cell's level-1 ancestor within its top-level face cell.
	 *
	 * @param int $level
	 * @return int
	 */
	public function childPosition(int $level): int {
		return (int)S2::unsignedRightShift($this->id, (2 * (self::MAX_LEVEL - $level) + 1)) & 3;
	}

	// Methods that return the range of cell ids that are contained
	// within this cell (including itself). The range is *inclusive*
	// (i.e. test using >= and <=) and the return values of both
	// methods are valid leaf cell ids.
	//
	// These methods should not be used for iteration. If you want to
	// iterate through all the leaf cells, call child_begin(MAX_LEVEL) and
	// child_end(MAX_LEVEL) instead.
	//
	// It would in fact be error-prone to define a range_end() method,
	// because (range_max().id() + 1) is not always a valid cell id, and the
	// iterator would need to be tested using "<" rather that the usual "!=".
	/**
	 * @return S2CellId
	 */
	public function rangeMin(): S2CellId {
		return new S2CellId($this->id - ($this->lowestOnBit() - 1));
	}

	/**
	 * @return S2CellId
	 */
	public function rangeMax(): S2CellId {
		return new S2CellId($this->id + ($this->lowestOnBit() - 1));
	}

	/**
	 * Return true if the given cell is contained within this one.
	 *
	 * @param S2CellId $other
	 * @return bool
	 */
	public function contains(S2CellId $other): bool {
		// assert (isValid() && other.isValid());
		return $other->greaterOrEquals($this->rangeMin()) && $other->lessOrEquals($this->rangeMax());
	}

	/**
	 * Return true if the given cell intersects this one.
	 *
	 * @param S2CellId $other
	 * @return bool
	 */
	public function intersects(S2CellId $other): bool {
		// assert (isValid() && other.isValid());
		return $other->rangeMin()->lessOrEquals($this->rangeMax()) && $other->rangeMax()->greaterOrEquals($this->rangeMin());
	}

	/**
	 * @param $level
	 * @return S2CellId
	 */
	public function parent($level = null): S2CellId {
		// assert (isValid() && level() > 0);
		if ($level === null) {
			$newLsb = $this->lowestOnBit() << 2;
		}
		else {
			$newLsb = self::lowestOnBitForLevel($level);
		}

		return new S2CellId(($this->id & -$newLsb) | $newLsb);
	}

	/**
	 * @param $level
	 * @return S2CellId
	 */
	public function childBegin($level = null): S2CellId {
		// assert (isValid() && level() < MAX_LEVEL);
		if ($level === null) {
			$oldLsb = $this->lowestOnBit();

			return new S2CellId($this->id - $oldLsb + ($oldLsb >> 2 & PHP_INT_MAX >> 1));
		}

		return new S2CellId($this->id - $this->lowestOnBit() + self::lowestOnBitForLevel($level));
	}

	/**
	 * @param $level
	 * @return S2CellId
	 */
	public function childEnd($level = null): S2CellId {
		if ($level === null) {
			// assert (isValid() && level >= this.level() && level <= MAX_LEVEL);
			return new S2CellId($this->id + $this->lowestOnBit() + self::lowestOnBitForLevel($level));
		}

		// assert (isValid() && level() < MAX_LEVEL);
		$oldLsb = $this->lowestOnBit();

		return new S2CellId($this->id + $oldLsb + S2::unsignedRightShift($oldLsb, 2));
	}

	// Iterator-style methods for traversing the immediate children of a cell or
	// all of the children at a given level (greater than or equal to the current
	// level). Note that the end value is exclusive, just like standard STL
	// iterators, and may not even be a valid cell id. You should iterate using
	// code like this:
	//
	// for(S2CellId c = id.childBegin(); !c.equals(id.childEnd()); c = c.next())
	// ...
	//
	// The convention for advancing the iterator is "c = c.next()", so be sure
	// to use 'equals()' in the loop guard, or compare 64-bit cell id's,
	// rather than "c !== id.childEnd()".

	/**
	 * Return the next cell at the same level along the Hilbert curve. Works
	 * correctly when advancing from one face to the next, but does *not* wrap
	 * around from the last face to the first or vice versa.
	 *
	 * @return S2CellId
	 */
	public function next(): S2CellId {
		return new S2CellId($this->id + ($this->lowestOnBit() << 1));
	}

	/**
	 * Return the previous cell at the same level along the Hilbert curve. Works
	 * correctly when advancing from one face to the next, but does *not* wrap
	 * around from the last face to the first or vice versa.
	 */
	public function prev(): S2CellId {
		return new S2CellId($this->id - ($this->lowestOnBit() << 1));
	}

	/**
	 * Like next(), but wraps around from the last face to the first and vice
	 * versa. Should *not* be used for iteration in conjunction with
	 * child_begin(), child_end(), Begin(), or End().
	 *
	 * @return S2CellId
	 */
	public function nextWrap(): S2CellId {
		$n = $this->next();
		if (self::unsignedLongLessThan($n->id, self::$WRAP_OFFSET)) {
			return $n;
		}

		return new S2CellId($n->id - self::$WRAP_OFFSET);
	}

	/**
	 * Like prev(), but wraps around from the last face to the first and vice
	 * versa. Should *not* be used for iteration in conjunction with
	 * child_begin(), child_end(), Begin(), or End().
	 *
	 * @return S2CellId
	 */
	public function prevWrap(): S2CellId {
		$p = $this->prev();
		if ($p->id < self::$WRAP_OFFSET) {
			return $p;
		}

		return new S2CellId($p->id + self::$WRAP_OFFSET);
	}

	/**
	 * @param int $level
	 * @return S2CellId
	 */
	public static function begin(int $level): S2CellId {
		return self::fromFacePosLevel(0, 0, 0)->childBegin($level);
	}

	/**
	 * @param int $level
	 * @return S2CellId
	 */
	public static function end(int $level): S2CellId {
		return self::fromFacePosLevel(5, 0, 0)->childEnd($level);
	}

	/**
	 * Decodes the cell id from a compact text string suitable for display or
	 * indexing. Cells at lower levels (i.e. larger cells) are encoded into
	 * fewer characters. The maximum token length is 16.
	 *
	 * @param string|null $token the token to decode
	 * @return S2CellId for that token
	 */
	public static function fromToken(string $token = null): S2CellId {
		if ($token === null) {
			throw new RuntimeException('Null string in S2CellId.fromToken');
		}
		if ($token === '') {
			throw new RuntimeException('Empty string in S2CellId.fromToken');
		}
		if ($token === 'X' || strlen($token) > 16) {
			return self::none();
		}

		// $value = hexdec(strrev($token));
		$value = hexdec($token);

		return new S2CellId($value);
	}

	/**
	 * Encodes the cell id to compact text strings suitable for display or indexing.
	 * Cells at lower levels (i.e. larger cells) are encoded into fewer characters.
	 * The maximum token length is 16.
	 *
	 * Simple implementation: convert the id to hex and strip trailing zeros. We
	 * could use base-32 or base-64, but assuming the cells used for indexing
	 * regions are at least 100 meters across (level 16 or less), the savings
	 * would be at most 3 bytes (9 bytes hex vs. 6 bytes base-64).
	 *
	 * @return string the encoded cell id
	 */
	public function toToken(): string {
		if ($this->id === 0) {
			return "X";
		}

		//$hex = Long.toHexString(id).toLowerCase(Locale.ENGLISH);
		$hex = strtolower(dechex($this->id));

		$sb = '';
		for ($i = strlen($hex); $i < 16; $i++) {
			$sb .= '0';
		}
		$sb .= $hex;

		for ($len = 16; $len > 0; $len--) {
			if ($sb[$len - 1] !== '0') {
				return substr($sb, 0, $len);
			}
		}

		throw new RuntimeException('Shouldn\'t make it here');
	}

	/**
	 * Returns true if (current * radix) + digit is a number too large to be
	 * represented by an unsigned long.  This is useful for detecting overflow
	 * while parsing a string representation of a number.
	 * Does not verify whether supplied radix is valid, passing an invalid radix
	 * will give undefined results or an ArrayIndexOutOfBoundsException.
	 *
	 * @param     $current
	 * @param int $digit
	 * @param int $radix
	 * @return bool
	 */
	private static function overflowInParse($current, int $digit, int $radix): bool {
		if ($current >= 0) {
			if ($current < self::$maxValueDivs[$radix]) {
				return false;
			}
			if ($current > self::$maxValueDivs[$radix]) {
				return true;
			}

			// current === maxValueDivs[radix]
			return ($digit > self::$maxValueMods[$radix]);
		}

		// current < 0: high bit is set
		return true;
	}

	// calculated as 0xffffffffffffffff / radix
	private static array $maxValueDivs = [
		0, 0, // 0 and 1 are invalid
		(float)'9223372036854775807L', (float)'6148914691236517205L', (float)'4611686018427387903L', // 2-4
		(float)'3689348814741910323L', (float)'3074457345618258602L', (float)'2635249153387078802L', // 5-7
		(float)'2305843009213693951L', (float)'2049638230412172401L', (float)'1844674407370955161L', // 8-10
		(float)'1676976733973595601L', (float)'1537228672809129301L', (float)'1418980313362273201L', // 11-13
		(float)'1317624576693539401L', (float)'1229782938247303441L', (float)'1152921504606846975L', // 14-16
		(float)'1085102592571150095L', (float)'1024819115206086200L', (float)'970881267037344821L', // 17-19
		(float)'922337203685477580L', (float)'878416384462359600L', (float)'838488366986797800L', // 20-22
		(float)'802032351030850070L', (float)'768614336404564650L', (float)'737869762948382064L', // 23-25
		(float)'709490156681136600L', (float)'683212743470724133L', (float)'658812288346769700L', // 26-28
		(float)'636094623231363848L', (float)'614891469123651720L', (float)'595056260442243600L', // 29-31
		(float)'576460752303423487L', (float)'558992244657865200L', (float)'542551296285575047L', // 32-34
		(float)'527049830677415760L', (float)'512409557603043100L' // 35-36
	];

	// calculated as 0xffffffffffffffff % radix
	private static array $maxValueMods = [
		0, 0, // 0 and 1 are invalid
		1, 0, 3, 0, 3, 1, 7, 6, 5, 4, 3, 2, 1, 0, 15, 0, 15, 16, 15, 15, // 2-21
		15, 5, 15, 15, 15, 24, 15, 23, 15, 15, 31, 15, 17, 15, 15 // 22-36
	];

	/**
	 * Return the four cells that are adjacent across the cell's four edges.
	 * Neighbors are returned in the order defined by S2Cell::GetEdge. All
	 * neighbors are guaranteed to be distinct.
	 *
	 * @param S2CellId $neighbors
	 */
	public function getEdgeNeighbors(S2CellId $neighbors): void {
		$i = 0;
		$j = 0;

		$level = $this->level();
		$size = 1 << (self::MAX_LEVEL - $level);
		$face = $this->toFaceIJOrientation($i, $j);

		// Edges 0, 1, 2, 3 are in the S, E, N, W directions.
		$neighbors[0] = self::fromFaceIJSame($face, $i, $j - $size, $j - $size >= 0)->parent($level);
		$neighbors[1] = self::fromFaceIJSame($face, $i + $size, $j, $i + $size < self::MAX_SIZE)->parent($level);
		$neighbors[2] = self::fromFaceIJSame($face, $i, $j + $size, $j + $size < self::MAX_SIZE)->parent($level);
		$neighbors[3] = self::fromFaceIJSame($face, $i - $size, $j, $i - $size >= 0)->parent($level);
	}

	/**
	 * Return the neighbors of closest vertex to this cell at the given level, by
	 * appending them to "output". Normally there are four neighbors, but the
	 * closest vertex may only have three neighbors if it is one of the 8 cube
	 * vertices.
	 *
	 * Requires: level < this.evel(), so that we can determine which vertex is
	 * closest (in particular, level === MAX_LEVEL is not allowed).
	 *
	 * @param $level
	 * @param $output
	 */
	public function getVertexNeighbors($level, &$output): void {
		// "level" must be strictly less than this cell's level so that we can
		// determine which vertex this cell is closest to.
		// assert (level < this.level());
		$i = 0;
		$j = 0;
		$face = $this->toFaceIJOrientation($i, $j);

		// Determine the i- and j-offsets to the closest neighboring cell in each
		// direction. This involves looking at the next bit of "i" and "j" to
		// determine which quadrant of this->parent(level) this cell lies in.
		$halfSize = 1 << (self::MAX_LEVEL - ($level + 1));
		$size = $halfSize << 1;
		if (($i & $halfSize) !== 0) {
			$iOffset = $size;
			$iSame = ($i + $size) < self::MAX_SIZE;
		}
		else {
			$iOffset = -$size;
			$iSame = ($i - $size) >= 0;
		}
		if (($j & $halfSize) !== 0) {
			$jOffset = $size;
			$jSame = ($j + $size) < self::MAX_SIZE;
		}
		else {
			$jOffset = -$size;
			$jSame = ($j - $size) >= 0;
		}

		$output[] = $this->parent($level);
		$output[] = self::fromFaceIJSame($face, $i + $iOffset, $j, $iSame)->parent($level);
		$output[] = self::fromFaceIJSame($face, $i, $j + $jOffset, $jSame)->parent($level);
		// If i- and j- edge neighbors are *both* on a different face, then this
		// vertex only has three neighbors (it is one of the 8 cube vertices).
		if ($iSame || $jSame) {
			$output[] = self::fromFaceIJSame($face, $i + $iOffset, $j + $jOffset, $iSame && $jSame)->parent($level);
		}
	}

	/**
	 * Return a leaf cell given its cube face (range 0..5) and i- and
	 * j-coordinates (see s2.h).
	 *
	 * @param $face
	 * @param $i
	 * @param $j
	 * @return S2CellId
	 */
	public static function fromFaceIJ($face, $i, $j): S2CellId {
		// Optimization notes:
		// - Non-overlapping bit fields can be combined with either "+" or "|".
		// Generally "+" seems to produce better code, but not always.

		// gcc doesn't have very good code generation for 64-bit operations.
		// We optimize this by computing the result as two 32-bit integers
		// and combining them at the end. Declaring the result as an array
		// rather than local variables helps the compiler to do a better job
		// of register allocation as well. Note that the two 32-bits halves
		// get shifted one bit to the left when they are combined.
		$n = array(0, $face << (self::POS_BITS - 33));

		// Alternating faces have opposite Hilbert curve orientations; this
		// is necessary in order for all faces to have a right-handed
		// coordinate system.
		$bits = ($face & self::SWAP_MASK);

		// Each iteration maps 4 bits of "i" and "j" into 8 bits of the Hilbert
		// curve position. The lookup table transforms a 10-bit key of the form
		// "iiiijjjjoo" to a 10-bit value of the form "ppppppppoo", where the
		// letters [ijpo] denote bits of "i", "j", Hilbert curve position, and
		// Hilbert curve orientation respectively.

		for ($k = 7; $k >= 0; --$k) {
			$bits = self::getBits($n, $i, $j, $k, $bits);
		}

		return new S2CellId(((($n[1] << 32) + $n[0]) << 1) + 1);
	}

	/**
	 * @param $n
	 * @param $i
	 * @param $j
	 * @param $k
	 * @param $bits
	 * @return int|mixed
	 */
	private static function getBits(&$n, $i, $j, $k, $bits) {
		$mask = (1 << self::LOOKUP_BITS) - 1;
		$bits += ((($i >> ($k * self::LOOKUP_BITS)) & $mask) << (self::LOOKUP_BITS + 2));
		$bits += ((($j >> ($k * self::LOOKUP_BITS)) & $mask) << 2);
		$bits = self::$LOOKUP_POS[$bits];
		$n[$k >> 2] |= (($bits >> 2) << (($k & 3) * 2 * self::LOOKUP_BITS));
		$bits &= (self::SWAP_MASK | self::INVERT_MASK);

		return $bits;
	}

	/**
	 * Return the (face, i, j) coordinates for the leaf cell corresponding to this
	 * cell id. Since cells are represented by the Hilbert curve position at the
	 * center of the cell, the returned (i,j) for non-leaf cells will be a leaf
	 * cell adjacent to the cell center. If "orientation" is non-NULL, also return
	 * the Hilbert curve orientation for the current cell.
	 *
	 * @param      $pi
	 * @param      $pj
	 * @param null $orientation
	 * @return int
	 */
	public function toFaceIJOrientation(&$pi, &$pj, &$orientation = null): int {
		// System.out.println("Entering toFaceIjorientation");
		$face = $this->face();
		$bits = ($face & self::SWAP_MASK);

		// System.out.println("face = " + face + " bits = " + bits);

		// Each iteration maps 8 bits of the Hilbert curve position into
		// 4 bits of "i" and "j". The lookup table transforms a key of the
		// form "ppppppppoo" to a value of the form "iiiijjjjoo", where the
		// letters [ijpo] represents bits of "i", "j", the Hilbert curve
		// position, and the Hilbert curve orientation respectively.
		//
		// On the first iteration we need to be careful to clear out the bits
		// representing the cube face.
		for ($k = 7; $k >= 0; --$k) {
			$bits = $this->getBits1($pi, $pj, $k, $bits);
			// System.out.println("pi = " + pi + " pj= " + pj + " bits = " + bits);
		}

		if ($orientation !== null) {
			// The position of a non-leaf cell at level "n" consists of a prefix of
			// 2*n bits that identifies the cell, followed by a suffix of
			// 2*(MAX_LEVEL-n)+1 bits of the form 10*. If n==MAX_LEVEL, the suffix is
			// just "1" and has no effect. Otherwise, it consists of "10", followed
			// by (MAX_LEVEL-n-1) repetitions of "00", followed by "0". The "10" has
			// no effect, while each occurrence of "00" has the effect of reversing
			// the kSwapMask bit.
			// assert (S2.POS_TO_ORIENTATION[2] === 0);
			// assert (S2.POS_TO_ORIENTATION[0] === S2.SWAP_MASK);
			if (($this->lowestOnBit() & 0x1111111111111110) !== 0) {
				$bits ^= self::SWAP_MASK;
			}
			$orientation = $bits;
		}

		return $face;
	}

	/**
	 * @param $i
	 * @param $j
	 * @param $k
	 * @param $bits
	 * @return int|mixed
	 */
	private function getBits1(&$i, &$j, $k, $bits) {
		$nbits = ($k === 7) ? (self::MAX_LEVEL - 7 * self::LOOKUP_BITS) : self::LOOKUP_BITS;

		$shift = ($k * 2 * self::LOOKUP_BITS + 1);

		$bits += (($this->id >> $shift & PHP_INT_MAX >> ($shift - 1)) & ((1 << (2 * $nbits)) - 1)) << 2;
		/* >>> */
		/*
		 * System.out.println("id is: " + id_); System.out.println("bits is " +
		 * bits); System.out.println("lookup_ij[bits] is " + lookup_ij[bits]);
		 */
		$bits = self::$LOOKUP_IJ[$bits];
		$i += ($bits >> (self::LOOKUP_BITS + 2)) << ($k * self::LOOKUP_BITS);
		/*
		* System.out.println("left is " + ((bits >> 2) & ((1 << kLookupBits) -
		* 1))); System.out.println("right is " + (k * kLookupBits));
		* System.out.println("j is: " + j.intValue()); System.out.println("addition
		* is: " + ((((bits >> 2) & ((1 << kLookupBits) - 1))) << (k *
		* kLookupBits)));
		*/
		$j += ((($bits >> 2) & ((1 << self::LOOKUP_BITS) - 1))) << ($k * self::LOOKUP_BITS);
		$bits &= (self::SWAP_MASK | self::INVERT_MASK);

		return $bits;
	}

	/** Return the lowest-numbered bit that is on for cells at the given level.
	 *
	 * @return int
	 */
	public function lowestOnBit(): int {
		return $this->id & -$this->id;
	}

	/**
	 * Return the lowest-numbered bit that is on for this cell id, which is equal
	 * to (uint64(1) << (2 * (MAX_LEVEL - level))). So for example, a.lsb() <=
	 * b.lsb() if and only if a.level() >= b.level(), but the first test is more
	 * efficient.
	 *
	 * @param $level
	 * @return int
	 */
	public static function lowestOnBitForLevel($level): int {
		return 1 << (2 * (self::MAX_LEVEL - $level));
	}

	/**
	 * Return the i- or j-index of the leaf cell containing the given s- or t-value.
	 *
	 * @param $s
	 * @return mixed
	 */
	private static function stToIJ($s) {
		// Converting from floating-point to integers via static_cast is very slow
		// on Intel processors because it requires changing the rounding mode.
		// Rounding to the nearest integer using FastIntRound() is much faster.

		$m = self::MAX_SIZE / 2; // scaling multiplier

		return max(0, min(2 * $m - 1, round($m * $s + ($m - 0.5))));
	}

	/**
	 * Convert (face, si, ti) coordinates (see s2.h) to a direction vector (not
	 * necessarily unit length).
	 *
	 * @param $face
	 * @param $si
	 * @param $ti
	 * @return S2Point
	 */
	private static function faceSiTiToXYZ($face, $si, $ti): S2Point {
		$kScale = 1.0 / self::MAX_SIZE;
		$u = S2Projections::stToUV($kScale * $si);
		$v = S2Projections::stToUV($kScale * $ti);

		return S2Projections::faceUvToXyz($face, $u, $v);
	}

	/**
	 * Given (i, j) coordinates that may be out of bounds, normalize them by
	 * returning the corresponding neighbor cell on an adjacent face.
	 *
	 * @param $face
	 * @param $i
	 * @param $j
	 * @return S2CellId
	 */
	private static function fromFaceIJWrap($face, $i, $j): S2CellId {
		// Convert i and j to the coordinates of a leaf cell just beyond the
		// boundary of this face. This prevents 32-bit overflow in the case
		// of finding the neighbors of a face cell, and also means that we
		// don't need to worry about the distinction between (s,t) and (u,v).
		$i = max(-1, min(self::MAX_SIZE, $i));
		$j = max(-1, min(self::MAX_SIZE, $j));

		// Find the (s,t) coordinates corresponding to (i,j). At least one
		// of these coordinates will be just outside the range [0, 1].
		$kScale = 1.0 / self::MAX_SIZE;
		$s = $kScale * (($i << 1) + 1 - self::MAX_SIZE);
		$t = $kScale * (($j << 1) + 1 - self::MAX_SIZE);

		// Find the leaf cell coordinates on the adjacent face, and convert
		// them to a cell id at the appropriate level.
		$p = S2Projections::faceUvToXyz($face, $s, $t);
		$face = S2Projections::xyzToFace($p);
		$st = S2Projections::validFaceXyzToUv($face, $p);

		return self::fromFaceIJ($face, self::stToIJ($st->getX()), self::stToIJ($st->getY()));
	}

	/**
	 * Public helper function that calls FromFaceIJ if sameFace is true, or
	 * FromFaceIJWrap if sameFace is false.
	 *
	 * @param $face
	 * @param $i
	 * @param $j
	 * @param $sameFace
	 * @return S2CellId
	 */
	public static function fromFaceIJSame($face, $i, $j, $sameFace): S2CellId {
		if ($sameFace) {
			return self::fromFaceIJ($face, $i, $j);
		}

		return self::fromFaceIJWrap($face, $i, $j);
	}

	/**
	 * @param $that
	 * @return bool
	 */
	public function equals($that): bool {
		if (!($that instanceof self)) {
			return false;
		}

		return $this->id() === $that->id();
	}

	/**
	 * Returns true if x1 < x2, when both values are treated as unsigned.
	 *
	 * @param $x1
	 * @param $x2
	 * @return bool
	 */
	public static function unsignedLongLessThan($x1, $x2): bool {
		return ($x1 + PHP_INT_MIN) < ($x2 + PHP_INT_MIN);
	}

	/**
	 * Returns true if x1 > x2, when both values are treated as unsigned.
	 *
	 * @param $x1
	 * @param $x2
	 * @return bool
	 */
	public static function unsignedLongGreaterThan($x1, $x2): bool {
		return ($x1 & ~PHP_INT_MAX) > ($x2 & ~PHP_INT_MAX);
	}

	/**
	 * @param S2CellId $x
	 * @return bool
	 */
	public function lessThan(S2CellId $x): bool {
		return self::unsignedLongLessThan($this->id, $x->id);
	}

	/**
	 * @param S2CellId $x
	 * @return bool
	 */
	public function greaterThan(S2CellId $x): bool {
		return self::unsignedLongGreaterThan($this->id, $x->id);
	}

	/**
	 * @param S2CellId $x
	 * @return bool
	 */
	public function lessOrEquals(S2CellId $x): bool {
		return self::unsignedLongLessThan($this->id, $x->id) || $this->id === $x->id;
	}

	/**
	 * @param S2CellId $x
	 * @return bool
	 */
	public function greaterOrEquals(S2CellId $x): bool {
		return self::unsignedLongGreaterThan($this->id, $x->id) || $this->id === $x->id;
	}

	/**
	 * @return int
	 * @Override
	 */
	public function hashCode(): int {
		return (int)(S2::unsignedRightShift($this->id, 32) + $this->id);
	}

	/**
	 * @return string
	 */
	public function __toString() {
		return sprintf('(face=%d, pos=%16x, level=%d)', $this->face(), $this->pos(), $this->level());
	}

	/**
	 * @param $level
	 * @param $i
	 * @param $j
	 * @param $origOrientation
	 * @param $pos
	 * @param $orientation
	 * @return void
	 */
	public static function initLookupCell($level, $i, $j, $origOrientation, $pos, $orientation): void {
		if ($level === self::LOOKUP_BITS) {
			$ij = ($i << self::LOOKUP_BITS) + $j;
			self::$LOOKUP_POS[($ij << 2) + $origOrientation] = ($pos << 2) + $orientation;
			self::$LOOKUP_IJ[($pos << 2) + $origOrientation] = ($ij << 2) + $orientation;
		}
		else {
			$level++;
			$i <<= 1;
			$j <<= 1;
			$pos <<= 2;
			// Initialize each sub-cell recursively.
			for ($subPos = 0; $subPos < 4; $subPos++) {
				$ij = S2::posToIJ($orientation, $subPos);
				$orientationMask = S2::posToOrientation($subPos);
				self::initLookupCell($level, $i + ($ij >> 1 & PHP_INT_MAX >> 0), $j + ($ij & 1), $origOrientation, $pos + $subPos, $orientation ^ $orientationMask);
				/* >>> */
			}
		}
	}

	/**
	 * @param S2CellId $that
	 * @return int
	 */
	public function compareTo(S2CellId $that): int {
		if (self::unsignedLongLessThan($this->id, $that->id)) {
			return -1;
		}

		return self::unsignedLongGreaterThan($this->id, $that->id) ? 1 : 0;
	}
}

S2CellId::$LOOKUP_POS = array_pad(array(), 1 << (2 * S2CellId::LOOKUP_BITS + 2), 0);
S2CellId::$LOOKUP_IJ = array_pad(array(), 1 << (2 * S2CellId::LOOKUP_BITS + 2), 0);
S2CellId::initLookupCell(0, 0, 0, 0, 0, 0);
S2CellId::initLookupCell(0, 0, 0, S2CellId::SWAP_MASK, 0, S2CellId::SWAP_MASK);
S2CellId::initLookupCell(0, 0, 0, S2CellId::INVERT_MASK, 0, S2CellId::INVERT_MASK);
S2CellId::initLookupCell(0, 0, 0, S2CellId::SWAP_MASK | S2CellId::INVERT_MASK, 0, S2CellId::SWAP_MASK | S2CellId::INVERT_MASK);
