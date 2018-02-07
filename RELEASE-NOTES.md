## ChangeLog v0.6.8-pre1 (20180207)

Minor release with a much faster binary. 10-20% faster than v0.6.6,
due to ldc and LLVM improvements. Fixes speed regression of v0.6.7 for
large files. See also [performance](https://github.com/biod/sambamba/blob/master/test/benchmark/stats.org)

+ Fixed Makefile for general use, see #332
+ Started benchmarking, see #283 and https://github.com/biod/sambamba/blob/master/test/benchmark/stats.org
+ Readded [Travis-ci support](https://travis-ci.org/biod/sambamba) for Linux (MacOS is disabled #338)
+ Documentation updates
