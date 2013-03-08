/++
	Merge Sort for Random-Access Ranges
	
	Written and tested for DMD 2.059 and Phobos
	
	Authors:  Xinok
	License:  Public Domain
	
	Bugs:
	Parallel sort fails to compile in debug builds

    -------------------------------
    March, 2013: added support for custom task pool //lomereiter
++/

module thirdparty.mergesort;
import std.range, std.algorithm, std.functional, std.array, std.parallelism;

/++
	Performs a merge sort on a random-access range according to predicate less.
	
	Returns: Sorted input as SortedRange
	
	Params:
	half = Set to true to merge using O(n/2) additional space, or false for O(n)
	task_pool = Thread pool to use
	temp = Optionally provide your own additional space for sorting
		
	Examples:
	-----------------
	int[] array = [10, 37, 74, 99, 86, 28, 17, 39, 18, 38, 70];
	mergeSort(array);
	mergeSort!"a > b"(array); // Sorts array descending	
    TaskPool pool = new TaskPool();
    scope(exit) pool.finish();
	mergeSort(array, pool);   // Sorts array using custom task pool
	
	int[] temp;
	temp.length = 64;
	mergeSort(array, pool, temp); // Sorts array using temporary memory provided by user
	-----------------
++/

@trusted SortedRange!(R, less) mergeSort(alias less = "a < b", bool half = true, R)(R range, TaskPool task_pool=taskPool, ElementType!(R)[] temp = null)
{
	static assert(isRandomAccessRange!R);
	static assert(hasLength!R);
	static assert(hasSlicing!R);
	static assert(hasAssignableElements!R);
	
	MergeSortImpl!(less, half, R).sort(range, task_pool, temp);
	
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
	void sort(R range, TaskPool task_pool, T[] temp)
	{
		static if(half)
		{
			if(temp.length < range.length / 2) temp.length = range.length / 2;
		}
		else
		{
			if(temp.length < range.length) temp.length = range.length;
		}
		
        concSort(range, task_pool, task_pool.size, temp);
	}
	
	/// Concurrently sort range
	void concSort(R range, TaskPool task_pool, size_t threadCount, T[] temp)
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
			auto th = task!(concSort)(range[0 .. mid], task_pool, threadCount / 2, temp[0 .. $ / 2]);
			task_pool.put(th);
			concSort(range[mid .. range.length], task_pool, threadCount - (threadCount / 2), temp[$ / 2 .. $]);
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
