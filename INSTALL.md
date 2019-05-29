# INSTALL SAMBAMBA

## Sambamba dependencies

* D compiler (ldc)
* BioD (git submodule)
* gcc tool chain (for htslib and lz4)
* htslib (git submodule)
* undeaD (git submodule)
* libz
* liblz4

## Install Sambamba from source

After checking out the source from github with git submodules is is
possibleto install the build tools with GNU Guix

    guix package -i gcc-toolchain gdb bash ld-wrapper ldc which python2 git

Even better, with Guix, you can create a light-weight container in the source tree
and run our development setup (gold was added lately by ldc)

    guix environment -C guix --ad-hoc gcc-toolchain gdb bash ld-wrapper ldc which python git binutils-gold vim
    make clean
    make -j 4
    make check

this way all dependencies are isolated.
