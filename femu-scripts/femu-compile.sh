#!/bin/bash

NRCPUS="$(cat /proc/cpuinfo | grep "vendor_id" | wc -l)"

make clean
# --disable-werror --extra-cflags=-w
../configure --extra-cflags=-Wnoimplicit-fallthrough --enable-kvm --target-list=x86_64-softmmu --with-git-submodules=validate
make -j $NRCPUS

echo ""
echo "===> FEMU compilation done ..."
echo ""
exit
