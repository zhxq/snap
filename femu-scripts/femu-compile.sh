#!/bin/bash

NRCPUS="$(cat /proc/cpuinfo | grep "vendor_id" | wc -l)"

make clean
# --disable-werror --extra-cflags=-w
# FEMU has some problem in their code - they did not add explicit fallthrough.
# We have to disable to default sanity check enabled in QEMU.
# There is no other way to override this in their configure script.
sed 's/-Wimplicit-fallthrough=2/-Wimplicit-fallthrough=0/g' ../configure > ../configure2
chmod +x ../configure2
../configure2 --enable-kvm --target-list=x86_64-softmmu --with-git-submodules=validate
rm ../configure2
make -j $NRCPUS

echo ""
echo "===> FEMU compilation done ..."
echo ""
exit
