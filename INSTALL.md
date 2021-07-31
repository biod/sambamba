# INSTALL SAMBAMBA

## Sambamba dependencies

* D compiler (ldc)
* gcc tool chain (for lz4)
* libz

## Install Sambamba from source

After checking out the source from github with git submodules it is
possible to install the build tools with GNU Guix

    guix package -i gcc-toolchain gdb bash ld-wrapper ldc which python git

Even better, with Guix, you can create a light-weight container in the source tree
and run our development setup (gold was added lately by ldc)

    guix environment -C guix --ad-hoc gcc-toolchain gdb bash ld-wrapper ldc which python git binutils-gold vim
    make clean
    env CC=gcc make -f Makefile.guix -j 4
    make -f Makefile.guix check

this way all dependencies are isolated. To create a static release use

    make static

## Development

We use GNU Guix containers for development. Install Guix and run a build
container with

    . .guix-build
    make -f Makefile.guix clean
    # build the debug version
    env CC=gcc make -f Makefile.guix lz4-static -j 8
    env CC=gcc make -f Makefile.guix -j 8
    make -f Makefile.guix check

To make the static release:

    env CC=gcc make -f Makefile.guix static

It gives some errors, but should work:

    ./bin/sambamba

When you only get unit tests disable them with `--DRT-testmode=run-main`

Note that this all also works in the emacs shell.

### Guix VM

    guix package -i qemu -p ~/opt/qemu
    . ~/opt/qemu/etc/profile

Download the bootable image from https://guix.gnu.org/en/download/ and
start it with, for example

    qemu-system-x86_64    -nic user,model=virtio-net-pci    -enable-kvm -m 1024    -device virtio-blk,drive=myhd    -drive if=none,file=guix-system-vm-image-1.2.0.x86_64-linux,id=myhd
