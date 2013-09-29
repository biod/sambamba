/++
	Unstable Sort for Random-Access Ranges
	
	Written and tested for DMD 2.058 and Phobos
	
	Authors:  Xinok
	License:  Public Domain
++/

module thirdparty.unstablesort;
import std.range, std.algorithm, std.functional, std.parallelism;

/++
	Performs an unstable sort on a random-access range according to predicate less.
	The algorithm is a quick sort which resorts to heap sort to avoid worst-case.
	
	Returns: Sorted input as SortedRange
	
	Params:
	range = Range to be sorted
	pool = Task pool to use
	
	Params:
	less = Predicate (string, function, or delegate) used for comparing elements; Defaults to "a < b"
	R = Type of range to be sorted; Must be a finite random-access range with slicing
	
	Examples:
	-----------------
	int[] array = [10, 37, 74, 99, 86, 28, 17, 39, 18, 38, 70];
	unstableSort(array, taskPool);   // Sorts array using multiple threads
	-----------------
++/

@trusted SortedRange!(R, less) unstableSort(alias less = "a < b", R)(R range, TaskPool pool)
{
	static assert(isRandomAccessRange!R);
	static assert(hasLength!R);
	static assert(hasSlicing!R);
	static assert(hasAssignableElements!R);
	
	UnstableSortImpl!(less, R).sort(range, pool);
	
	assert(isSorted!(less)(range.save), "Range is not sorted");
	return assumeSorted!(less, R)(range.save);
}

/// Unstable sort implementation
template UnstableSortImpl(alias pred, R)
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
	
	/// Entry sort function
	void sort(R range, TaskPool pool)
	{
		concSort(range, range.length, pool);
	}
	
	/// Recursively partition list
	void sort(R range, real depth, TaskPool pool)
	{
		while(true)
		{
			if(range.length <= MAX_INSERT)
			{
				binaryInsertionSort(range);
				return;
			}
			if(depth < 1.0)
			{
				heapSort(range);
				return;
			}
			
			depth /= 1.5;
			
			immutable mid = partition(range);
			
			if(mid <= range.length / 2)
			{
				sort(range[0 .. mid - 1], depth, pool);
				range = range[mid .. range.length];
			}
			else
			{
				sort(range[mid .. range.length], depth, pool);
				range = range[0 .. mid - 1];
			}
		}
	}
	
	/// Concurrently sorts range
	void concSort(R range, real depth, TaskPool pool)
	{
		if(range.length < MIN_THREAD)
		{
			sort(range, depth, pool);
			return;
		}
		if(depth < 1.0)
		{
			heapSort(range);
			return;
		}
		
		depth /= 1.5;
		
		immutable mid = partition(range);
		
		auto th = task!(concSort)(range[0 .. mid - 1], depth, pool);
		pool.put(th);
		concSort(range[mid .. range.length], depth, pool);
		th.workForce();
	}
	
	/// Partitions range, returns starting index of second range excluding pivot
	size_t partition(R range)
	{
		// Get median of five
		immutable b = range.length / 4, c = range.length / 2, d = b + c;
		medianSort(range[0], range[b], range[c], range[d], range[range.length - 1]);
		
		// Move first elements into place
		swap(range[1], range[b]);
		swap(range[2], range[c]);
		swap(range[range.length - 2], range[d]);
		
		// Variables
		T piv = range[2], o;
		size_t lef = 3, rig = range.length - 3;
		
		// Partition range
		while(lef < rig)
		{
			if(lessEqual(range[lef], piv)) ++lef;
			else
			{
				o = range[lef];
				range[lef] = range[rig];
				range[rig] = o;
				--rig;
			}
			
			// Checking for equality on both sides ensures a balanced
			// distribution of equal elements in both partitions
			if(greaterEqual(range[rig], piv)) --rig;
			else
			{
				o = range[lef];
				range[lef] = range[rig];
				range[rig] = o;
				++lef;
			}
		}
		
		 // This step is necessary to ensure pivot is inserted at correct location
		if(lessEqual(range[lef], piv)) ++lef;
		// Move pivot into place
		swap(range[lef - 1], range[2]);
		
		return lef;
	}
	
	/// Finds the median of five in six comparisons while satisfiying the condition: 
	/// (a < c && b < c && c < d && d < e)
	void medianSort(ref T a, ref T b, ref T c, ref T d, ref T e)
	out
	{
		assert(lessEqual(a, c) && lessEqual(b, c) && lessEqual(c, d) && lessEqual(c, e));
	}
	body
	{
		T o;
		
		if(greater(a, b)) swap(a, b);
		if(greater(d, e)) swap(d, e);
		
		if(greater(a, d))
		{
			o = a;
			a = d;
			d = b;
			b = o;
			if(greater(c, e))
			{
				o = c;
				c = d;
				d = e;
				e = o;
			}
			else swap(c, d);
		}
		else if(greater(b, c)) swap(b, c);
		
		if(greater(b, d))
		{
			swap(b, d);
			swap(c, e);
		}
		
		if(greater(c, d)) swap(c, d);
	}
	
	/// A simple insertion sort used for sorting small sublists
	void binaryInsertionSort(R range)
	{
		size_t lower, upper, center;
		T o;
		for(size_t i = 1; i < range.length; ++i)
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
	
	/// Bottom-up binary heap sort is used to avoid the worst-case of quick sort
	void heapSort(R range)
	{
		// Build Heap
		size_t i = (range.length - 2) / 2 + 1;
		while(i > 0) sift(range, --i, range.length);
		
		// Sort
		i = range.length - 1;
		while(i > 0)
		{
			swap(range[0], range[i]);
			sift(range, 0, i);
			--i;
		}
	}
	
	void sift(R range, size_t parent, immutable size_t end)
	{
		immutable root = parent;
		T value = range[parent];
		size_t child = void;
		
		// Sift down
		while(true)
		{
			child = parent * 2 + 1;
			
			if(child >= end) break;
			
			if(child + 1 < end && less(range[child], range[child + 1])) child += 1;
			
			range[parent] = range[child];
			parent = child;
		}
		
		child = parent;
		
		// Sift up
		while(child > root)
		{
			parent = (child - 1) / 2;
			if(less(range[parent], value))
			{
				range[child] = range[parent];
				child = parent;
			}
			else break;
		}
		
		range[child] = value;
	}
}

unittest
{
	bool testSort(alias pred, R)(R range)
	{
		unstableSort!(pred, R)(range, taskPool);
		return isSorted!pred(range);
	}
	
	int testCall(T)(in T[] arr)
	{
		int failures = 0;
		
		if(!testSort!"a < b"(arr.dup)) ++failures;
		if(!testSort!"a > b"(arr.dup)) ++failures;
		
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
}
