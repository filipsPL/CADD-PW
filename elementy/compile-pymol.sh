#!/bin/bash

cd pymol

prefix=/opt/pymol-svn
modules=$prefix/modules

# enable c++11
export CPPFLAGS="-std=c++0x"

# If you want to install as root, then split this line up in "build"
# and "install" and run the "install" with "sudo"
sudo python setup.py build install \
    --home=$prefix \
    --install-lib=$modules \
    --install-scripts=$prefix
