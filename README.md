[![Anaconda-Server Badge](https://anaconda.org/bioconda/sambamba/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.13200.svg)](http://dx.doi.org/10.5281/zenodo.13200)
# Sambamba

Sambamba is a high performance modern robust and fast tool (and
library), written in the D programming language, for working with SAM
and BAM files.  Current parallelised functionality is an important
subset of samtools functionality, including view, index, sort,
markdup, and depth. 

Because of efficient use of modern multicore CPUs, usually `sambamba` is much faster
than `samtools`. For example, indexing a 2.5 Gb BAM file (fully cached into RAM) 
on a 8 core machine utilizes all cores at 64% CPU:

    time sambamba index merged_NIT20120138_F3_20130715.bam -t8

      real    0m17.398s
      user    1m25.841s
      sys     0m3.752s

meanwhile samtools is *4x* slower:

    time samtools index merged_NIT20120138_F3_20130715.bam

      real    1m8.083s
      user    1m6.640s
      sys     0m1.448s

In practice, the speedup is usually smaller since I/O becomes a bottleneck.
Even so, it makes a big difference, shifting the focus to I/O optimization, i.e.
less temporary files, more UNIX pipes, faster disk storage, tweaking filesystem, etc.
Most tools in `sambamba` support piping: just specify `/dev/stdin` or `/dev/stdout` as filenames.

Notice that `samtools` implements parallel BAM compression in `sort` and `merge`, 
but `sambamba` should be faster for these tasks (given same amount of memory) 
due to more cache-friendly approach to parallelization.
If it is not the case for you, please file a bug.

Sambamba is free and open source software, licensed under GPLv2+.
See manual pages [online](https://lomereiter.github.io/sambamba/docs/sambamba-view.html) 
to know more about what is available and how to use it.

For more information on Sambamba you can contact Artem Tarasov and Pjotr Prins.

# Binaries

If you are already a Conda user, please use [`bioconda`](https://bioconda.github.io/) channel.
If you are not, strongly consider becoming one.

Users of Homebrew can also use the formula from `homebrew-science`.

For those not in the mood to learn/install new package managers, there are of course Github [releases](https://github.com/lomereiter/sambamba/releases).

# Compiling Sambamba

The preferred method for compiling Sambamba is with the LDC compiler
which targets LLVM.

## Compiling for Linux

The LDC compiler's github repository also provides binary images. The current
preferred release for sambamba is LDC - the LLVM D compiler (>= 0.16.1). After
installing LDC:

```sh
    git clone --recursive https://github.com/lomereiter/sambamba.git
    cd sambamba
    make sambamba-ldmd2-64
```

Installing LDC only means unpacking an archive and setting some environmental variables, e.g. unpacking into `$HOME`:
```sh
cd
wget https://github.com/ldc-developers/ldc/releases/download/v0.17.1/ldc2-0.17.1-linux-x86_64.tar.xz
tar xJf ldc2-0.17.1-linux-x86_64.tar.xz
export PATH=~/ldc2-0.17.1-linux-x86_64/bin/:$PATH
export LIBRARY_PATH=~/ldc2-0.17.1-linux-x86_64/lib/
```

## Compiling for Mac OS X

```sh
    brew install ldc
    git clone --recursive https://github.com/lomereiter/sambamba.git
    cd sambamba
    make sambamba-ldmd2-64
```

# Development

Sambamba development and issue tracker is on
[github](https://github.com/lomereiter/sambamba). Developer
documentation can be found in the source code and the [development
documentation](https://github.com/lomereiter/sambamba-dev-docs).

# Copyright

Sambamba is distributed under GNU Public License v2+.

# Citation

If you are using Sambamba in your research, please cite the following article:

A. Tarasov, A. J. Vilella, E. Cuppen, I. J. Nijman, and P. Prins. Sambamba: fast processing of NGS alignment formats. Bioinformatics, 2015.
