## ChangeLog v0.6.8 (20181003)

Pre-release with a much faster statically compiled binary. 10-20%
faster than v0.6.6, due to ldc and LLVM improvements. Fixes speed
regression of v0.6.7 for large files due to singleobj compilation. See
also #345 and
[performance](https://github.com/biod/sambamba/blob/master/test/benchmark/stats.org)

64-bit compilation should be fine on ldc 1.10. i386 target is still a problem.

+ Fix mark duplicates in files with many contigs, see #361 (thanks Devon Ryan @dpryan79)
+ Fix missing PM tag in #356 (thanks Kurt Hetrick @Kurt-Hetrick)
+ Fix Bcftools version checking #352 (thanks Nathan S. Watson-Haigh @nathanhaigh)
+ Fixate version info in BAM output headers for reproducibility. See #357
+ Fixed Makefile for general use, see #332
+ Started benchmarking, see #283 and https://github.com/biod/sambamba/blob/master/test/benchmark/stats.org
+ Readded [Travis-ci support](https://travis-ci.org/biod/sambamba) for Linux (MacOS is disabled #338)
+ Fixed MacOS build in Travis with ae269cfbdf2e78750ce7f8dc70ad32f80a6682df
+ Updated BioD to latest https://github.com/biod/BioD/commit/5e56b2bb45324af2194b3339d298fd827c8003ae
+ Bug fixes:
  * #328 Debug version: SAM output of CRAM file is populated with debug on pipe
  * #331 Segmentation fault attempting to view header in json format
  * #335 Intel Xeon bug may segfault Sambamba - this was tracked down to an Intel Xeon bug
  * #345 sambamba index 0.6.7 takes 4x longer than 0.6.6 on the same files
+ Documentation updates
+ Updated lz4 to latest
+ Added support for GNU Guix and build containers
+ Added shunit2 to the source tree for testing
+ Update python build dependencies to use python3.x
+ Fixed a number of D compiler messages on deprecated features (ldc 1.11)
