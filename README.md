The project aims to develop robust and fast library for working with BAM files.

# Tools

There're already some tools developed on top of the library, that are [faster](https://github.com/lomereiter/sambamba/wiki/Comparison-with-samtools) than samtools
in many cases. See `CLItools/` directory and manual pages on [wiki][] to see what is available and how to use it.

# Library features

### Parallelism
	
BGZF blocks are unpacked in parallel. 
This allows to loop over alignments 2-5 times faster, depending on how many cores you have.

Also, it's easy to implement parallel algorithms in D with <code>std.parallelism</code> module. 
But be aware that currently D garbage collector is of the stop-the-world type,
so you should avoid allocating memory on the heap in parallel loops. 

### You don't pay for what you don't use

Alignments are parsed lazily. So lazily that they're not parsed at all 
at the moment of their creation from the chunk of memory. 
That's possible because BAM is a binary format.

All access is done via _properties_ which quickly calculate where 
the requested field lies in the chunk, and return it to you. That works
not only for primitive types but for strings and arrays as well, thanks to
D _slices_ making it extremely easy.

### Native speed

The speed of D code compiled with GDC compiler is about the same as speed of C code, because
gcc and gdc share the same optimization backend. 
The code compiled with DMD runs slower, sometimes by a factor of 2, but the
advantages you get with it are more cool language features 
(like uniform function call syntax introduced in version 2.059) and faster compiling times. 

That being said, DMD is recommended for development, and GDC for faster execution times.
It is GDC which is used for building Debian packages.

### Small memory footprint

Most methods return ranges, not arrays. 

Ranges in D resemble iterators in Java, or yield statement in Python and Ruby.
They can be easily combined in various ways, allowing you to write code in functional style.

Typically, the compiler will do a lot of inlining work, so that the performance will be
about the same as of manually written deeply nested loops, while the code will remain clean.

### Tight integration with D language constructs

As said above, ranges are used almost everywhere. 
That means all loops can be written using <code>foreach</code> statement.

But that's not the only syntax sugar you can use. A simple example:

    import bamfile;
    import std.stdio;

	void main() {
        foreach(alignment; BamFile("mybamfile.bam")["chr1"][10_000 .. 12_000]) {
            write(alignment.read_name); 
            Value read_group = alignment["RG"];
            writeln(read_group.is_nothing ? "" : " " ~ to!string(read_group));
        }
    }

This code opens BAM file, fetches alignments with reference 'chr1' overlapping
[10_000, 12_000) interval (in zero-based coordinate system), and
prints read name and contents of 'RG' tag if available. Easy, huh? :-)

<code>Value</code> is a special type for representing tag values which is internally
just a tagged union, thus the overhead is minimal. 

### Clean and documented code with lots of unittests.

Unittests are yet another D programming language feature. 
You can be sure that the behaviour of the library won't change over the time,
except some bugs might go away ;-).

----------------------------------------------------------------------------------------------

##You're welcome to contribute and promote the library!

----------------------------------------------------------------------------------------------

#How to use

1. Clone the repository.
2. Download DMD compiler from http://dlang.org/download.html
3. Install the [Ragel][] state machine compiler.
4. Read GitHub [wiki][] of the project.

[wiki]: https://github.com/lomereiter/sambamba/wiki
[Ragel]: http://www.complang.org/ragel/

# Copyright

This library and tools that come with it are all distributed under GNU Public License v2+.

