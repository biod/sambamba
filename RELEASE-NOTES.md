## ChangeLog v0.8.2 (2022....)

+ Purged lz4 sources if favour of using the lz4 dependency (the Debian
  library appears to be complete
+ Fixed meson builds for Debian, see https://github.com/biod/sambamba/issues/484
  (thanks @atille)
+ Removed most ldc2 deprecation warnings
+ Moved Python scripts 2->3 (following the Debian patch)
+ Fix ldc 1.27.1 build by using get to unwrap Nullable
+ Merged PR #480 (thanks @dukc)

## ChangeLog v0.8.1 (20210731)

+ Some FASTA unittests break when they write to disk and read after -
  marked with ~version(Broken)~
+ Fix build for LDC 1.26.0 and LLVM 9.0.1 - some speed improvement

## ChangeLog v0.8.0 (20201130)

Maintainance release and bug fixes: this is a special release where we
removed all CRAM support. The added value of CRAM in sambamba was
limited because it was using essentially the same htslib backend as
samtools. Removing the htslib dependency removes one maintenance
headache. See also https://github.com/biod/sambamba/issues/425.

BioD was also moved back into the main trunk. We separated it in the
past, but as there is no development there we might as well have it in
Sambamba (again).

+ Removed CRAM support and htslib dependency with ec78eabfbdebd620cf5d4a4d8232692184eacbf7
+ Removed BioD source tree back into Sambamba
+ Sambamba builds on ARM64 (see [Travis-CI](https://travis-ci.org/github/biod/sambamba))
+ Changed github issue tracker template to report bugs only
+ Added FASTA handlers, indexing and slicing (thanks @NickRoz1)
+ Updated slice.d (thanks @NickRoz1)
+ Fixed FASTA file handlers not closing cdc5f84c5b6135eedeaf5cba265e3eb9a3c23227
+ Added badges to the top of README
+ Added Meson build system for Debian (thanks Matthias @ximion)
+ Fixed Python2 -> Python3 stuff (thanks Matthias @ximion)

## ChangeLog v0.7.1 (20191128)

+ Fixed major segfault problem that got introduced after ldc 1.10, see https://github.com/biod/sambamba/issues/393
+ make error message more informative - see https://github.com/biod/sambamba/pull/411/files, thanks https://github.com/EwaMarek

## ChangeLog v0.7.0 (20190529)

+ Release avoids crashing problem by building with ldc-1.10.0, see https://github.com/biod/sambamba/issues/393
+ Merge pull request #396 from joelmartin/flagstat-return returning 0

## ChangeLog v0.6.9 (20190313)

+ Added support for Picard style sorting, see https://github.com/biod/sambamba/issues/369 - thanks https://github.com/TimurIs
+ Add a new sorting option to pull together mates when sorting by read name, see https://github.com/biod/sambamba/pull/380 - thanks https://github.com/emi80
+ Many module renames to build a BioD for the future - thanks https://github.com/george-githinji

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
ee61000bcb33a82013c284bac8feb91f  sambamba-0.6.8.gz
gzip -d sambamba-0.6.8.gz
chmod a+x sambamba-0.6.8

./sambamba-0.8.6

sambamba 0.6.8 by Artem Tarasov and Pjotr Prins (C) 2012-2018
    LDC 1.10.0 / DMD v2.080.1 / LLVM6.0.1 / bootstrap LDC - the LLVM D compiler (0.17.4)
```

The binary images were built on x86_64 with

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
make clean && make static -j 16 && make check
```

Git submodule versions were

```
 git submodule status
 2f0634b187e0f454809432093238cf31e9fbfee6 BioD (v0.2.0-5-g2f0634b)
 2f3c3ea7b301f9b45737a793c0b2dcf0240e5ee5 htslib (0.2.0-rc10-271-g2f3c3ea)
 b3692db46d2b23a7c0af2d5e69988c94f126e10a lz4 (v1.8.2)
 9be93876982b5f14fcca60832563b3cd767dd84d undeaD (v1.0.1-49-g9be9387)
 ```
