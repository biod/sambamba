# Sambamba

Sambamba is a high performance modern robust and fast tool (and
library), written in the D programming language, for working
with BAM files.  Current functionality is an important subset of
samtools functionality. 

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

# Quick install from source

Download LDC2 compiler from Github: https://github.com/ldc-developers/ldc/releases
Add it to `$PATH`, fetch sambamba sources and run `make`, e.g. for Linux 64-bit:

    wget https://github.com/ldc-developers/ldc/releases/download/v0.12.0/ldc2-0.12.0-linux-x86_64.tar.xz
    tar xJf ldc2-0.12.0-linux-x86_64.tar.xz
    git clone --recursive https://github.com/lomereiter/sambamba.git
    cd sambamba
    export PATH=../ldc2-0.12.0-linux-x86_64/bin:$PATH
    make sambamba-ldmd2-64

The binary will be in `build/` subdirectory.

# Binaries

See Github [releases](https://github.com/lomereiter/sambamba/releases)

## Note

If you are going to build LDC compiler from source, add `-O3` flag in 
`build/runtime/CMakeFiles/phobos-ldc.dir/flags.make` when building LDC, 
otherwise Zlib library will be compiled without optimizations.

# Copyright

Sambamba is distributed under GNU Public License v2+.
