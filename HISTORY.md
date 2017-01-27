Sambamba does not see much development nowadays, and here I'll attempt to explain why, 
putting aside standard excuses such as 'lack of time' or 'lack of funding' but rather cutting straight to the _technical_ core of the matter.

## Original motivation

When this project was started (2012), there were few options to work with raw sequencing data. 
The choice was essentially between samtools and Picard, and both shared the huge drawback: they didn't scale.
They were written with an old single-threaded mindset, which is quite hard to get rid of.

Obviously, this was a very real pain point to address, and it turned out that not much is needed to have an edge: 
BAM files are composed of little compressed blocks, and (de-)compressing those in a multithreaded fashion is 
not that difficult in any _modern_ language. (note: C doesn't count as such, and neither does Java before version 8)
Back then we chose D for the job, and it wasn't a terrible choice, although the language still doesn't enjoy much popularity.

Even in 2012 there was at least one library that was well optimized for multithreaded BAM reading/writing, namely
[SeqAn](https://github.com/seqan/seqan), a very solid C++ library.
It was little known, however—presumably due to lack of useful user-facing command line tools such as those that come with samtools.
So having a set of operations reasonably close to the widely used samtools was the second necessary ingredient for a success.

## BAM format drawbacks

Although quite good, the format is far from perfect. 
First of all, the compression algorithm, gzip, is set in stone.
The original concept of samtools includes piping commands as much as possible, using no compression for intermediate steps.
This alleviates the problem somewhat, but there are lots of interesting alternatives, 
providing for example fast compression at the expense of some disk space 
(in particular [Snappy](https://github.com/google/snappy) and [lz4](https://github.com/lz4/lz4/), both appeared in 2011).

Another issue is that most tools access only a few fields, but nevertheless have to decompress whole record contents. 
That again leads to lots of computational power wasted on unnecessary decompression.

[CRAM](http://www.ebi.ac.uk/ena/software/cram-toolkit) introduced both columnar storage and compression algorithm 
flexibility, solving both aforementioned issues. Its developers attempt to _"carefully balance support for utility 
(e.g streaming, indexing, direct computational access) with sufficiently deep compression"_, and with reference compression 
it is indeed a great choice for long-term storage.

However, in 2013 happened something game-changing: two flexible and efficient columnar storage formats, 
[Parquet](https://parquet.apache.org/) and [Orc](https://orc.apache.org/), appeared on the stage and quickly gained 
in popularity afterwards, becoming standard in big data processing communities.

And when we introduce sustainability into the CRAM balance equation, it might well be that scarce on human resources 
computational genomics community would be better off by adopting one of these widely used formats, so that it can enjoy
library and tooling support coming for free from the thriving big data industry.

## Wrong bet on language

Before Sambamba came into being, Pjotr Prins had written an interesting [article](http://thebird.nl/blog/D_Dragon.html) 
trying to figure out what is the best programming language for bioinformatics. He was clearly struggling to make 
the choice between Scala and D, weighed many factors in and came to the conclusion that D is more performant because "it is much easier to get low level tweaked code".

I would like to argue that at least _since 2015_ it is much easier to get _high level_ yet still _tweaked_ code in *Scala*. 
What does it even mean, you ask? Go and read [this article on Project Tungsten](https://databricks.com/blog/2015/04/28/project-tungsten-bringing-spark-closer-to-bare-metal.html), which freed [Spark](https://spark.apache.org/) of its main performance issues: the new architecture avoids JVM object overhead and generates efficient cache-aware code on the fly. Moreover, IBM engineers are [working](https://developer.ibm.com/javasdk/2016/12/22/how-the-jit-compiler-can-exploit-simd-and-gpu-for-spark-workloads/) on improving the code generator to use SIMD instructions. Now realize that a typical bioinformatician never went beyond Bash and Python (ok, maybe Perl), and doesn't have any clue as to the meaning of 'SIMD', 'JVM object overhead' or 'cache-aware', so you certainly cannot expect/trust them to do any low-level tweaking!

Today Scala (+Spark) lets bioinformaticians focus on biology and algorithms, providing out of the box decent performance 
and opportunity to parallelize data processing on the cluster. (Although Python users can use PySpark, it suffers from 
overhead of JVM<->Python object conversion, which will hopefully be mitigated by [Apache Arrow](http://arrow.apache.org/) 
in the next few years, so stay tuned and learn some Scala meanwhile.)

## ADAM

All provided links should have convinced you that Parquet/ORC are wonderful storage formats, and Scala is a great language. 
Combine the puzzle, and what you get is [ADAM](https://github.com/bigdatagenomics/adam), building upon Spark and Parquet 
and following best practices of software development. It wouldn't have been possible 5 years ago, without 
the strong foundations provided by open-source community, but now that these are available it's highly beneficial 
to go ahead and embrace them. I must point out, however, that it was only in 2015 that Spark performance was given 
a drastic boost, so it must have been a pretty shaky foundation on the start, but luckily turned out quite well.

Just a couple of years ago I was surprised to discover that ADAM authors called Sambamba a 'legacy tool' in 
a [paper](https://dl.acm.org/citation.cfm?id=2742787). Since then I came to realize why it is indeed legacy, 
and this document serves to shed some light on this fact. ADAM is an example of solid engineering, and it also takes 
scaling to the whole next level, i.e. beyond a single node. So try it out, check if it works for you, and fill their
bug tracker with your use cases as their team is heading towards 1.0 release in the next year.

&nbsp;

Artem Tarasov, January 2017
