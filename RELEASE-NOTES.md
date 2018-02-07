## ChangeLog v0.6.8-pre1 (20180207)

Minor release with a much faster binary. 10-20% faster than v0.6.6,
due to ldc and LLVM improvements. Fixes speed regression of v0.6.7 for
large files. See also [performance](https://github.com/biod/sambamba/blob/master/test/benchmark/stats.org)

+ Fixed Makefile for general use, see #332
+ Started benchmarking, see #283 and https://github.com/biod/sambamba/blob/master/test/benchmark/stats.org
+ Readded [Travis-ci support](https://travis-ci.org/biod/sambamba) for Linux (MacOS is disabled #338)
+ Updated BioD to latest https://github.com/biod/BioD/commit/5e56b2bb45324af2194b3339d298fd827c8003ae
+ Bug fixes:
  * #328 Debug version: SAM output of CRAM file is populated with debug on pipe
  * #331 Segmentation fault attempting to view header in json format
  * #335 Intel Xeon bug may segfault Sambamba - this was tracked down to an Intel Xeon bug
+ Documentation updates
