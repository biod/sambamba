## ChangeLog v0.6.8 (20181004)

Pre-release with a much faster statically compiled binary. 10-20%
faster than v0.6.6, due to ldc and LLVM improvements. Fixes speed
regression of v0.6.7 for large files due to singleobj compilation. See
also #345 and
[performance](https://github.com/biod/sambamba/blob/master/test/benchmark/stats.org)

64-bit compilation should be fine on ldc 1.10+. i386 target is still a problem.

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
+ Updated lz4 to latest (still source in tree because Debian dropped frame support in liblz4-dev)
+ Added support for GNU Guix and build containers
+ Added shunit2 to the source tree for testing
+ Update python build dependencies to use python3.x
+ Fixed a number of D compiler messages on deprecated features (ldc 1.11)

To install the image, download and

```sh
md5sum sambamba-0.6.8.gz
25efb5604ae5fe7c750e8020326787c5  sambamba-0.8.6.gz
gzip -d sambamba-0.6.8.gz
chmod a+x sambamba-0.6.8

./sambamba-0.8.6

sambamba 0.6.8 by Artem Tarasov and Pjotr Prins (C) 2012-2018
    LDC 1.10.0 / DMD v2.080.1 / LLVM6.0.1 / bootstrap LDC - the LLVM D compiler (0.17.4)
```

The binary images were reproducibly built on x86_64 with

```sh
~/.config/guix/current/bin/guix pull -l
Generation 3    Sep 25 2018 09:39:08
  guix 932839f
    repository URL: https://git.savannah.gnu.org/git/guix.git
    branch: origin/master
    commit: 932839ff124ff3b0dd3070914fb1c5beec69bf32


guix environment -C guix --ad-hoc gcc gdb bash ld-wrapper ldc which python git
make clean && make -j 16 && make check

for x in `ldd bin/sambamba|cut -d ' ' -f 3` ; do realpath $x ; done
/gnu/store/l4lr0f5cjd0nbsaaf8b5dmcw1a1yypr3-glibc-2.27/lib/libpthread-2.27.so
/gnu/store/l4lr0f5cjd0nbsaaf8b5dmcw1a1yypr3-glibc-2.27/lib/libm-2.27.so
/gnu/store/l4lr0f5cjd0nbsaaf8b5dmcw1a1yypr3-glibc-2.27/lib/librt-2.27.so
/gnu/store/l4lr0f5cjd0nbsaaf8b5dmcw1a1yypr3-glibc-2.27/lib/libdl-2.27.so
/gnu/store/bmaxmigwnlbdpls20px2ipq1fll36ncd-gcc-8.2.0-lib/lib/libgcc_s.so.1
/gnu/store/l4lr0f5cjd0nbsaaf8b5dmcw1a1yypr3-glibc-2.27/lib/libc-2.27.so
# build static image
make clean && make release -j 16 && make check
```

Git submodule versions were

```
 git submodule status
 2f0634b187e0f454809432093238cf31e9fbfee6 BioD (v0.2.0-5-g2f0634b)
 2f3c3ea7b301f9b45737a793c0b2dcf0240e5ee5 htslib (0.2.0-rc10-271-g2f3c3ea)
 b3692db46d2b23a7c0af2d5e69988c94f126e10a lz4 (v1.8.2)
 9be93876982b5f14fcca60832563b3cd767dd84d undeaD (v1.0.1-49-g9be9387)
 ```
