#!/bin/bash
# Huaicheng Li <huaicheng@cs.uchicago.edu>
# Run FEMU as a black-box SSD (FTL managed by the device)

# image directory
IMGDIR=$HOME/images
# Virtual machine disk image
OSIMGF=$IMGDIR/u20s.qcow2

if [[ ! -e "$OSIMGF" ]]; then
	echo ""
	echo "VM disk image couldn't be found ..."
	echo "Please prepare a usable VM image and place it as $OSIMGF"
	echo "Once VM disk image is ready, please rerun this script again"
	echo ""
	exit
fi

sudo x86_64-softmmu/qemu-system-x86_64 \
    -name "FEMU-BBSSD-VM" \
    -enable-kvm \
    -cpu host \
    -smp 8 \
    -m 24G \
    -device virtio-scsi-pci,id=scsi0 \
    -device scsi-hd,drive=hd0 \
    -drive file=$OSIMGF,if=none,aio=native,cache=none,format=qcow2,id=hd0 \
    -device femu,devsz_mb=18432,femu_mode=1,stream_support=on,max_streams=16,death_time_prediction=on \
    -net user,hostfwd=tcp::8080-:22,hostfwd=tcp::18081-:8081,hostfwd=tcp::18082-:8082,hostfwd=tcp::18083-:8083,hostfwd=tcp::18084-:8084 \
    -net nic,model=virtio \
    -nographic \
    -virtfs local,path=$(mkdir -p ../hostshare && cd ../hostshare && pwd),mount_tag=host0,security_model=passthrough,id=host0 \
    -device virtio-9p-pci,fsdev=host0,mount_tag=hostshare  \
    -qmp unix:./qmp-sock,server,nowait 2>&1 | tee log
