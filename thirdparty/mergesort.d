/++
	Merge Sort for Random-Access Ranges
	
	Written and tested for DMD 2.059 and Phobos
	
	Authors:  Xinok
	License:  Public Domain
	
	Bugs:
	Parallel sort fails to compile in debug builds
++/

module thirdparty.mergesort;
import std.range, std.algorithm, std.functional, std.array, std.parallelism;

/++
	Performs a merge sort on a random-access range according to predicate less.
	
	Returns: Sorted input as SortedRange
	
	Params:
	half = Set to true to merge using O(n/2) additional space, or false for O(n)
	threaded = Set to true to sort using multiple threads
	temp = Optionally provide your own additional space for sorting
		
	Examples:
	-----------------
	int[] array = [10, 37, 74, 99, 86, 28, 17, 39, 18, 38, 70];
	mergeSort(array);
	mergeSort!"a > b"(array); // Sorts array descending	
	mergeSort(array, true);   // Sorts array using multiple threads
	
	int[] temp;
	temp.length = 64;
	mergeSort(array, false, temp); // Sorts array using temporary memory provided by user
	-----------------
++/

@trusted SortedRange!(R, less) mergeSort(alias less = "a < b", bool half = true, R)(R range, bool threaded = false, ElementType!(R)[] temp = null)
{
	static assert(isRandomAccessRange!R);
	static assert(hasLength!R);
	static assert(hasSlicing!R);
	static assert(hasAssignableElements!R);
	
	MergeSortImpl!(less, half, R).sort(range, threaded, temp);
	
	if(!__ctfe) assert(isSorted!(less)(range.save), "Range is not sorted");
	return assumeSorted!(less, R)(range.save);
}

/// Merge Sort implementation
template MergeSortImpl(alias pred, bool half, R)
{
	static assert(isRandomAccessRange!R);
	static assert(hasLength!R);
	static assert(hasSlicing!R);
	static assert(hasAssignableElements!R);
		
	alias ElementType!R T;
	
	alias binaryFun!pred less;
	bool greater(T a, T b){ return less(b, a); }
	bool greaterEqual(T a, T b){ return !less(a, b); }
	bool lessEqual(T a, T b){ return !less(b, a); }
	
	enum MAX_INSERT = 32;        // Maximum length for an insertion sort
	enum MIN_THREAD = 1024 * 64; // Minimum length of a sublist to initiate new thread

	/// Entry point for merge sort
	void sort(R range, bool threaded, T[] temp)
	{
		static if(half)
		{
			if(temp.length < range.length / 2) temp.length = range.length / 2;
		}
		else
		{
			if(temp.length < range.length) temp.length = range.length;
		}
		
		if(threaded && !__ctfe)
			concSort(range, defaultPoolThreads + 1, temp);
		else
			split(range, temp);
	}
	
	/// Concurrently sort range
	void concSort(R range, size_t threadCount, T[] temp)
	{
		if(threadCount < 2 || range.length < MIN_THREAD)
		{
			split(range, temp);
			return;
		}
		
		debug
		{
			//@ Threading code currently does not compile in debug builds
			split(range, temp);
		}
		else
		{
			immutable mid = range.length / 2;
			auto th = task!(concSort)(range[0 .. mid], threadCount / 2, temp[0 .. $ / 2]);
			taskPool.put(th);
			concSort(range[mid .. range.length], threadCount - (threadCount / 2), temp[$ / 2 .. $]);
			th.workForce();
			merge(range, mid, temp);
		}
	}

	
	/// Recursively split range and merge halves
	void split(R range, T[] temp)
	{
		if(range.length <= MAX_INSERT)
		{
			binaryInsertionSort(range);
			return;
		}
		immutable mid = range.length / 2;
		split(range[0 .. mid], temp);
		split(range[mid .. range.length], temp);
		merge(range, mid, temp);
	}
		
	/// Merge two halves using temp
	static if(half)
	void merge(R range, immutable size_t mid, T[] temp)
	{
		assert(mid <= range.length);
		assert(temp.length >= range.length / 2);
		
		temp = temp[0 .. mid];
		copy(range[0..mid], temp);
		
		size_t i = 0, lef = 0, rig = mid;
		
		while(true)
		{
			if(lessEqual(temp[lef], range[rig]))
			{
				range[i++] = temp[lef++];
				if(lef >= temp.length) return;
			}
			else
			{
				range[i++] = range[rig++];
				if(rig >= range.length) while(true)
				{
					range[i++] = temp[lef++];
					if(lef >= temp.length) return;
				}
			}
		}
	}

	static if(!half)
	void merge(R range, immutable size_t mid, T[] temp)
	{
		assert(mid <= range.length);
		assert(temp.length >= range.length);
		
		size_t i = 0, lef = 0, rig = mid;
		while(true)
		{
			if(lessEqual(range[lef], range[rig]))
			{
				temp[i++] = range[lef++];
				if(lef >= mid) break;
			}
			else
			{
				temp[i++] = range[rig++];
				if(rig >= range.length)
				{
					while(lef < mid) temp[i++] = range[lef++];
					break;
				}
			}
		}
		copy(temp[0 .. i], range[0 .. i]);
	}
	
	/// A simple insertion sort used for sorting small sublists
	void binaryInsertionSort(R range)
	{
		size_t lower, upper, center;
		T o;
		foreach(i; 0 .. range.length)
		{
			o = range[i];
			lower = 0;
			upper = i;
			while(upper != lower)
			{
				center = (lower + upper) / 2;
				if(less(o, range[center])) upper = center;
				else lower = center + 1;
			}
			for(upper = i; upper > lower; --upper) range[upper] = range[upper-1];
			range[upper] = o;
		}
	}
	
	//@ Workaround for DMD issue 7898
	static if(__VERSION__ == 2059)
	void copy(R1, R2)(R1 src, R2 dst)
	{
		import std.traits;
		static if(isArray!R1 && isArray!R2) if(__ctfe)
		{
			dst[] = src[];
			return;
		}
		std.algorithm.copy(src, dst);
	}
}

unittest
{
	bool testSort(alias pred, bool half = false, R)(R range)
	{
		mergeSort!(pred, half, R)(range);
		return isSorted!pred(range);
	}
	
	int testCall(T)(in T[] arr)
	{
		int failures = 0;
		
		// Sort using O(n) space
		if(!testSort!("a < b", false)(arr.dup)) ++failures;
		if(!testSort!("a > b", false)(arr.dup)) ++failures;
		
		// Sort using O(n/2) space
		if(!testSort!("a < b", true)(arr.dup)) ++failures;
		if(!testSort!("a > b", true)(arr.dup)) ++failures;
		
		return failures;
	}
	
	// Array containing 256 random ints
	enum test = [
		10, 37, 74, 99, 86, 28, 17, 39, 18, 38, 70, 89, 94, 32, 46, 76, 43, 33, 62, 76, 
		37, 93, 45, 48, 49, 21, 67, 56, 58, 17, 15, 41, 91, 94, 95, 41, 38, 80, 37, 24, 
		26, 71, 87, 54, 72, 60, 29, 37, 41, 99, 31, 66, 75, 72, 86, 97, 37, 25, 98, 89, 
		53, 45, 52, 76, 51, 38, 59, 53, 74, 96, 94, 42, 68, 84, 65, 27, 49, 57, 53, 74, 
		39, 75, 39, 26, 46, 37, 68, 96, 19, 79, 73, 83, 36, 90, 11, 39, 48, 94, 97, 72, 
		37, 43, 69, 36, 41, 47, 31, 48, 33, 21, 20, 18, 45, 28, 47, 54, 41, 28, 47, 44, 
		51, 15, 21, 64, 82, 23, 41, 82, 30, 25, 78, 72, 50, 34, 45, 59, 14, 71, 50, 97, 
		39, 87, 74, 60, 52, 17, 87, 45, 69, 54, 91, 68, 46, 99, 78, 33, 27, 53, 41, 84, 
		82, 54, 29, 55, 53, 87, 13, 98, 55, 33, 73, 64, 19, 81, 57, 78, 23, 45, 94, 75, 
		55, 43, 93, 85, 96, 82, 44, 73, 22, 79, 89, 20, 36, 11, 12, 51, 86, 86, 75, 66, 
		81, 90, 80, 80, 36, 36, 47, 43, 86, 96, 45, 73, 70, 90, 57, 23, 86, 29, 12, 54, 
		37, 17, 87, 12, 36, 78, 26, 28, 30, 15, 10, 53, 76, 34, 23, 49, 65, 17, 37, 51, 
		26, 23, 66, 12, 26, 84, 60, 47, 30, 26, 78, 20, 42, 40, 63, 40
	];
	
	// Runtime test
	assert(testCall(test) == 0);
	
	// CTFE Test
	{
		enum result = testCall(test);
		static if(result != 0) pragma(msg, __FILE__, "(", __LINE__, "): Warning: mergeSort CTFE unittest failed ", result, " of 4 tests");
	}
	
	// Stability test
	bool icmp(ubyte a, ubyte b)
	{
		if(a >= 'a') a -= 'a' - 'A';
		if(b >= 'a') b -= 'a' - 'A';
		return a < b;
	}
	ubyte[] str = cast(ubyte[])"ksugnqtoyedwpvbmifaclrhjzxWELPGDVJIHBAMZCFUNORKSTYXQ".dup;
	mergeSort!icmp(str);
	assert(str == "aAbBcCdDeEfFgGhHiIjJkKlLmMnNoOpPqQrRsStTuUvVwWxXyYzZ");
}
