#!/bin/bash
# Huaicheng <huaicheng@cs.uchicago.edu>
# Please run this script as root.

SYSTEM=`uname -s`

if [[ -f /etc/debian_version ]]; then
	# Includes Ubuntu, Debian
    apt-get install -y gcc pkg-config git libglib2.0-dev libfdt-dev libpixman-1-dev zlib1g-dev
    apt-get install -y libaio-dev

	# Additional dependencies
	apt-get install -y libnuma-dev
    apt-get install -y ninja-build
    apt-get install -y libcap-ng-dev libattr1 libattr1-dev
else
    echo "pkgdep: unsupported system type ($SYSTEM), please install QEMU depencies manually"
	exit 1
fi

mkdir -p ../hostshare

echo "===> Dependency installation ... Done!"
echo "===> You can put things you want to share between host and virtual machine"
echo "======> in the \"hostshare\" folder in your FEMU code base ($(cd ../hostshare && pwd))."