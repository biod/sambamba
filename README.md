# Sambamba

Sambamba is a high performance modern robust and fast tool (and
library), written in the D programming language, for working
with BAM files.  Current functionality is an important subset of
samtools functionality. 

Because of efficient use of modern multicore CPUs, usually Sambamba is much faster
than samtools. For example, indexing an 18 Gb BAM file on a fast 8 core
machine utilizes all cores at 45% CPU:

    time ~/sambamba index /scratch/HG00119.mapped.ILLUMINA.bwa.GBR.exome.20111114.bam            
      real    1m42.930s
      user    6m19.964s
      sys     0m32.362s

meanwhile samtools

    time ~/samtools index /scratch/HG00119.mapped.ILLUMINA.bwa.GBR.exome.20111114.bam 
      real    5m37.669s
      user    5m9.127s
      sys     0m13.605s

Such a speedup can make a difference dealing with 1000 genomes (or more).

See also a further
[comparison](https://github.com/lomereiter/sambamba/wiki/Comparison-with-samtools)
on more limited hardware.

Sambamba is free and open source software. Sambamba filtering is part of the Galaxy [tool
shed](http://toolshed.g2.bx.psu.edu/repos/lomereiter/sambamba_filter)
and will be part of [CloudBiolinux](http://cloudbiolinux.org/) soon.
A Debian package is available for download.

See manual pages on [wiki](https://github.com/lomereiter/sambamba/wiki) to know more about 
what is available and how to use it.

For more information on Sambamba you can contact Artem Tarasov and Pjotr Prins.

# Quick install

The only dependency is a D2 compiler (dmd >= 2.062). E.g. from the source
tree on 64-bit Debian/Ubuntu

    git clone --recursive https://github.com/lomereiter/sambamba
    wget http://ftp.digitalmars.com/dmd_2.062-0_amd64.deb
    sudo dpkg -i dmd_2.062-0_amd64.deb
    make

where the binaries reside in ./build

# Binaries for 64-bit systems

https://www.dropbox.com/sh/v05fsb5aarob3xe/iUHgyud31a/sambamba

## Note

Since version 0.3.0, release builds are done using [LDC2 compiler](http://github.com/ldc-developers/ldc)

If you are going to build executable yourself, add `-O3` flag in 
`build/runtime/CMakeFiles/phobos-ldc.dir/flags.make` when building LDC, 
otherwise Zlib library will be compiled without optimizations.

# Copyright

Sambamba is distributed under GNU Public License v2+.
