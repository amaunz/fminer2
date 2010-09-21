#!/bin/bash

./liblast/configure
./libbbrc/configure
./fminer/configure

make -C liblast/
make -C libbbrc/
make -C fminer/
