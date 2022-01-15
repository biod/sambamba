# INSTALL SAMBAMBA

## Sambamba dependencies

* D compiler (ldc)
* libz
* liblz4
* meson + ninja for builds

The biod project was moved back into the sambamba git repo. It is no longer a dependency. We use [.guix-build](.guixbuild) to set up the build environment in a GNU Guix container.

## Install Sambamba from source

Tested on Debian:

    meson build --buildtype release
    cd build
    ninja
    ninja test

After checking out the source from github with git submodules it is
possible to install the build tools with GNU Guix

    guix package -i gcc-toolchain gdb bash ld-wrapper ldc which python git

Even better, with Guix, you can create a light-weight container in the source tree
and run our development setup (gold was added lately by ldc)

    source .guix-build
    make clean
    env CC=gcc make -f Makefile.guix -j 4
    make -f Makefile.guix check ./bin/sambamba-0.8.2

this way all dependencies are isolated. To create a static release use

    env CC=gcc make -f Makefile.guix static -j 4 (FIXME)

Alternatively use the meson+ninja build with

    rm -rf build/ ; env D_LD=gold CC=gcc meson build --buildtype release
    cd build/
    env CC=gcc ninja
    env CC=gcc ninja test

## Development

We use GNU Guix containers for development. Install Guix and run a build
container with

    . .guix-build
    make -f Makefile.guix clean
    # Set versions
    python3 ./gen_ldc_version_info.py ldc2 > utils/ldc_version_info_.d
    # build the debug version
    env CC=gcc make -f Makefile.guix lz4-static -j 8
    env CC=gcc make -f Makefile.guix -j 8

Instead, to make the release:

    env CC=gcc make -f Makefile.guix -j 8 release

To make the static release:

    env CC=gcc make -f Makefile.guix static

Run tests

    make -f Makefile.guix check

Run binary

    ./bin/sambamba

When you only get unit tests disable them with `--DRT-testmode=run-main`

Note that this all also works in the emacs shell.

### Guix VM

    guix package -i qemu -p ~/opt/qemu
    . ~/opt/qemu/etc/profile

Download the bootable image from https://guix.gnu.org/en/download/ and
start it with, for example

    qemu-system-x86_64    -nic user,model=virtio-net-pci    -enable-kvm -m 1024    -device virtio-blk,drive=myhd    -drive if=none,file=guix-system-vm-image-1.2.0.x86_64-linux,id=myhd

## Check list for release

- [ ] Test meson build with local lz4
- [ ] Build and test static version
- [ ] Build and test optimized version with performance run
- [ ] Update release notes
- [ ] Release on github
- [ ] Ping Debian, Guix and Conda projects
