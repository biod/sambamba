[![Build Status](https://travis-ci.org/biod/sambamba.svg?branch=master)](https://travis-ci.org/biod/sambamba) [![DL](https://anaconda.org/bioconda/sambamba/badges/downloads.svg)](https://anaconda.org/bioconda/sambamba) [![BrewBadge](https://img.shields.io/badge/%F0%9F%8D%BAbrew-sambamba-brightgreen.svg)](https://github.com/brewsci/homebrew-bio) [![GuixBadge](https://img.shields.io/badge/gnuguix-sambamba-brightgreen.svg)](https://packages.guix.gnu.org/packages/sambamba/) [![DebianBadge](https://badges.debian.net/badges/debian/testing/sambamba/version.svg)](https://packages.debian.org/testing/sambamba)

A. Tarasov, A. J. Vilella, E. Cuppen, I. J. Nijman, and P. Prins. [Sambamba: fast processing of NGS alignment formats](https://doi.org/10.1093/bioinformatics/btv098). Bioinformatics, 2015.

# SAMBAMBA

Table of Contents
=================

   * [Introduction](#introduction)
   * [Binary Installation](#binary-installation)
      * [Install Stable Release](#install-stable-release)
      * [Bioconda Install](#bioconda-install)
      * [GNU Guix Install](#gnu-guix-install)
      * [Debian GNU/Linux Install](#debian-gnulinux-install)
      * [Homebrew Install](#homebrew-install)
   * [Getting Help](#getting-help)
      * [Reporting a Sambamba Bug or Issue](#reporting-a-sambamba-bug-or-issue)
      * [Checklist](#checklist)
      * [Code of Conduct](#code-of-conduct)
   * [Compiling Sambamba](#compiling-sambamba)
      * [Compilation Dependencies](#compilation-dependencies)
      * [Compiling for Linux](#compiling-for-linux)
         * [GNU Guix](#gnu-guix)
      * [Compiling for Mac OS X](#compiling-for-mac-os-x)
      * [Development](#development)
   * [Debugging and Troubleshooting](#debugging-and-troubleshooting)
      * [Segfaults on Certain Intel Xeons](#segfaults-on-certain-intel-xeons)
      * [Dump Core](#dump-core)
      * [Use catchsegv](#use-catchsegv)

* [Using GDB](#using-gdb)
   * [License](#license)
   * [Credit](#credit)


<a name="intro"></a>
# Introduction

*January 2023*: announcing v1.0 release of the great sambamba tool!

A minor fix and a major release. After 10 years and over one thousand citations we can announce sambamba 1.0 stable!

Sambamba is a high performance, highly parallel, robust and fast tool (and library), written in the D programming language, for working with SAM and BAM files. Because of its efficiency, Sambamba is an important work horse running in many sequencing centers around the world today. As of January 2023, Sambamba has been cited over [1000 times](http://scholar.google.nl/citations?hl=en&user=5ijHQRIAAAAJ) and has been installed from Conda over [250,000 times](https://anaconda.org/bioconda/sambamba). Sambamba is also distributed by [Debian](https://packages.debian.org/testing/sambamba). To cite Sambamba, see [Credit](#credit).

Current functionality is an important subset of `samtools`
functionality including view, index, sort, markdup, and depth. Most
tools support piping. Just specify `/dev/stdin` or `/dev/stdout` as
filenames. When we started writing Sambamba in 2012, the main
advantage over `samtools` was parallelized BAM reading and writing.

In March 2017, `samtools` 1.4 was released, reaching some parity. A [recent performance comparison](https://github.com/guigolab/sambamBench-nf) shows that Sambamba still holds its ground and can even do better. Here are some comparison [metrics](https://public-docs.crg.es/rguigo/Data/epalumbo/sambamba_ws_report.html). For example, for flagstat, Sambamba is 1.4x faster than samtools. For index, they are similar. For markdup, almost *6x* faster, and for view, *4x* faster. For sort, Sambamba has been beaten, though Sambamba is notably up to *2x* faster than samtools on large RAM machines (120GB+).

In addition, Sambamba has a few interesting features to offer. In particular:

* Fast large machine `sort`, see [performance](./test/benchmark/stats.org)
* Automatic index creation when writing any coordinate-sorted file
* `view -L <bed file>` utilizes BAM index to skip unrelated chunks
* `depth` allows to measure base, sliding window, or region coverages
  * [Chanjo](https://www.chanjo.co/) builds upon this and gets you to exon/gene levels of abstraction
* `markdup`, a fast implementation of Picard algorithm
* `slice` quickly extracts a region into a new file, tweaking only first/last chunks
* ...and more (you'll have to try)

The D language is extremely suitable for high performance computing (HPC). At this point, we think that the BAM format is here to stay for processing reference guided sequencing data, and we aim to make it easy to parse and process BAM files.

Sambamba is free and open source software, licensed under GPLv2+.
See manual pages [online](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
to know more about what is available and how to use it.

For more information on Sambamba, contact the mailing list (see [Getting help](#getting-help)).

## No CRAM Support

**Important Notice: With version 0.8 support for CRAM was removed from
Sambamba (see the [RELEASE NOTES](./RELEASE-NOTES.md))**

To use CRAM, you can still use one of the older (binary)
[releases](https://github.com/biod/sambamba/releases) of Sambamba.

<a name="install"></a>
# Binary Installation
## Install Stable Release

For those not in the mood to learn/install new package managers, there
are Github source and binary
[releases](https://github.com/biod/sambamba/releases). Simply download
the tarball, unpack it, and run it according to the accompanying
release notes.

Below package managers Conda, GNU Guix, Debian and Homebrew also
provide recent binary installs for Linux.  For MacOS you may use Conda
or Homebrew.

## Bioconda Install

[![Install with CONDA](https://anaconda.org/bioconda/sambamba/badges/installer/conda.svg)](https://anaconda.org/bioconda/sambamba)

There should be binary downloads for Linux and MacOS.

With Conda use the [`bioconda`](https://bioconda.github.io/) channel.

## GNU Guix Install

[![GuixBadge](https://img.shields.io/badge/gnuguix-sambamba--0.6.8-brightgreen.svg)](https://www.gnu.org/software/guix/packages/S/)

A reproducible [GNU Guix package](https://www.gnu.org/software/guix/packages/s.html) for Sambamba is available. The development version is packaged [here](https://gitlab.com/genenetwork/guix-bioinformatics/blob/master/gn/packages/sambamba.scm).

## Debian GNU/Linux Install

[![DebianBadge](https://badges.debian.net/badges/debian/testing/sambamba/version.svg)](https://packages.debian.org/testing/sambamba)

See also Debian package [status](https://tracker.debian.org/pkg/sambamba).

## Homebrew Install

[![BrewBadge](https://img.shields.io/badge/%F0%9F%8D%BAbrew-sambamba-brightgreen.svg)](https://github.com/brewsci/homebrew-bio)

Users of Homebrew can also use the formula from [homebrew-bio](https://github.com/brewsci/homebrew-bio).

    brew install brewsci/bio/sambamba

It should work for Linux and MacOS.

<a name="help"></a>
# Getting Help

Sambamba has a
[mailing list](https://groups.google.com/forum/#!forum/sambamba-discussion)
for installation help and general discussion.

## Reporting a Ssambamba Bug or Issue

Before posting an issue, search the issue tracker and [mailing list](https://groups.google.com/forum/#!forum/sambamba-discussion) first. It is likely someone may have encountered something similar. Also try running the latest version of Sambamba to make sure it has not been fixed already. Support/installation questions should be aimed at the mailing list. The issue tracker is for development issues around the software itself. When reporting an issue, include the output of the program and the contents of the output directory.

Please use the following checklist. It exists for multiple reasons. :)

## Checklist

1. [X] I have found an issue with Sambamba
2. [ ] I have searched for it on the [issue tracker](https://github.com/biod/sambamba/issues) (also check closed issues)
3. [ ] I have searched for it on the [mailing list](https://groups.google.com/forum/#!forum/sambamba-discussion)
4. [ ] I have tried the latest [release](https://github.com/biod/sambamba/releases) of Sambamba
5. [ ] I have read and agreed to the below code of conduct
6. [ ] If it is a support/install question, I have posted it to the [mailing list](https://groups.google.com/forum/#!forum/sambamba-discussion)
7. [ ] If it is software development related, I have posted a new issue on the [issue tracker](https://github.com/biod/sambamba/issues) or added to an existing one
8. [ ] In the message, I have included the output of my Sambamba run
9. [ ] In the message, I have included relevant files in the output directory
10. [ ] (Optional) I have made the data available to reproduce the problem

To find bugs, the Sambamba software developers may ask to install a development version of the software. They may also ask you for your data and will treat it confidentially.  Please always remember that Sambamba is written and maintained by volunteers with good intentions. Our time is valuable, too. By helping us as much as possible, we can provide this tool for everyone to use.

## Code of Conduct

By using Sambamba and communicating with its community, you implicitly agree to abide by the [code of conduct](https://software-carpentry.org/conduct/) as published by the Software Carpentry initiative.


<a name="compile"></a>
# Compiling Sambamba

Note: In general, there is no need to compile Sambamba. You can use a recent binary install as listed above.

The preferred method for compiling Sambamba is with the LDC compiler which targets LLVM. LLVM version 6 is faster than earlier editions.

## Compilation Dependencies

See [INSTALL.md](./INSTALL.md).

## Compiling for Linux

The LDC compiler's Github repository provides binary images. The current preferred release for Sambamba is LDC - the LLVM D compiler (>= 1.6.1). After installing LDC from https://github.com/ldc-developers/ldc/releases/ with, for example

```sh
cd
wget https://github.com/ldc-developers/ldc/releases/download/v$ver/ldc2-1.7.0-linux-x86_64.tar.xz
tar xvJf ldc2-1.7.0-linux-x86_64.tar.xz
export PATH=$HOME/ldc2-1.7.0-linux-x86_64/bin:$PATH
export LIBRARY_PATH=$HOME/ldc2-1.7.0-linux-x86_64/lib
```

```sh
git clone --recursive https://github.com/biod/sambamba.git
cd sambamba
make
```

To build a development/debug version run

```sh
make clean && make debug
```

To run the test, fetch shunit2 from https://github.com/kward/shunit2 and put it in the path so
you can run

```sh
make check
```

See also [INSTALL.md](./INSTALL.md).

### GNU Guix

Our development and release environment is GNU Guix.  To build
sambamba the LDC compiler is also available in GNU Guix:

```sh
guix package -A ldc
```

For more instructions see [INSTALL.md](./INSTALL.md).

## Compiling with Meson/Ninja

Debian uses a Meson+Ninja build. It may work with something like

```sh
meson build
cd build
ninja
ninja test

time ./build/sambamba sort HG00100.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam
sambamba 1.0.1
 by Artem Tarasov and Pjotr Prins (C) 2012-2023
    LDC 1.32.0 / DMD v2.102.2 / LLVM14.0.6 / bootstrap LDC - the LLVM D compiler (1.32.0)
real    0m13.343s
user    2m11.663s
sys     0m4.232s
```

or possibly with some additional tuning runtimes get close to the optimized static build (see [benchmarks](test/benchmark/stats.org).

```sh
rm -rf build/ ; env D_LD=gold CC=gcc meson build --buildtype release
cd build/
env CC=gcc ninja
env CC=gcc ninja test

time ./build/sambamba sort HG00100.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam
sambamba 1.0.1
 by Artem Tarasov and Pjotr Prins (C) 2012-2023
    LDC 1.32.0 / DMD v2.102.2 / LLVM14.0.6 / bootstrap LDC - the LLVM D compiler (1.32.0)
real    0m10.227s
user    2m7.203s
sys     0m4.039s
```

## Compiling for MacOS

Sambamba builds on MacOS. We have a Travis [integration test](https://travis-ci.org/pjotrp/sambamba) as
an example. It can be something like

```sh
    brew install ldc
    git clone --recursive https://github.com/biod/sambamba.git
    cd sambamba
    make
```

## Development

Sambamba development and issue tracker is on
[github](https://github.com/biod/sambamba). Developer
documentation can be found in the source code and the [development
documentation](https://github.com/biod/sambamba-dev-docs).

<a name="debug"></a>

# Debugging and Troubleshooting

## Segfaults on Certain Intel Xeons

Important Note: Some older Xeon processors segfault under heavy hyper
threading - which Sambamba utilizes.  Please read
[this](https://blog.cloudflare.com/however-improbable-the-story-of-a-processor-bug/)
when encountering seemingly random crashes. There is no real fix other
than disabling hyperthreading. Also discussed [here](https://github.com/biod/sambamba/issues/335). Thank Intel for producing this bug.

## Dump Core

In a crash, Sambamba can dump a core file. To make this happen set

```sh
ulimit -c unlimited
```

and run your command. Send us the core file so we can reproduce the state at
time of segfault.

## Use catchsegv

Another option is to use `catchsegv`

```sh
catchsegv ./build/sambamba command
```

This will show state on stdout which can be sent to us.

## Using GDB

In case of crashes, it's helpful to have GDB stacktraces (`bt`
command). A full stacktrace for all threads:

```
thread apply all backtrace full
```

Note that GDB should be made aware of D garbage collector, which emits
SIGUSR signals and GDB needs to ignore them with

```
handle SIGUSR1 SIGUSR2 nostop noprint
```

A binary, relocatable install of Sambamba with debug information and
all dependencies can be fetched from the binary link above.  Unpack
the tarball and run the contained install.sh script with TARGET

```
./install.sh ~/sambamba-test
```

Run Sambamba in GDB with

```
gdb -ex 'handle SIGUSR1 SIGUSR2 nostop noprint' \
  --args ~/sambamba-test/sambamba-*/bin/sambamba view --throw-error
```

<a name="license"></a>
# License

Sambamba is generously distributed under GNU Public License v2+.

<a name="credits"></a>
# Credit

Citations are the bread and butter of science.  If you are using Sambamba in your research and want to support our future work on Sambamba, please cite the following publication:
A. Tarasov, A. J. Vilella, E. Cuppen, I. J. Nijman, and P. Prins. [Sambamba: fast processing of NGS alignment formats](https://doi.org/10.1093/bioinformatics/btv098). Bioinformatics, 2015.

## Bibtex reference

```bibtex

@article{doi:10.1093/bioinformatics/btv098,
  author = {Tarasov, Artem and Vilella, Albert J. and Cuppen, Edwin and Nijman, Isaac J. and Prins, Pjotr},
  title = {Sambamba: fast processing of NGS alignment formats},
  journal = {Bioinformatics},
  volume = {31},
  number = {12},
  pages = {2032-2034},
  year = {2015},
  doi = {10.1093/bioinformatics/btv098},
  URL = { + http://dx.doi.org/10.1093/bioinformatics/btv098}
```
