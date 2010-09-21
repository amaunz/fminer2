#!/bin/bash
# Try running this file for a complete build

make -C liblast clean >/dev/null 2>&1
make -C libbbrc clean >/dev/null 2>&1
make -C fminer clean >/dev/null 2>&1

cd ./liblast/ ; ./configure >/dev/null
if [ $? -ne 0 ]; then
  echo "Configuration failed."
  exit 1
fi
cd -
make -C liblast/ > /dev/null 2>&1
if [ $? -ne 0 ]; then
  echo "Build failed."
  exit 1
fi


cd ./libbbrc/ ; ./configure >/dev/null
if [ $? -ne 0 ]; then
  echo "Configuration failed."
  exit 1
fi
cd -
make -C libbbrc/ > /dev/null 2>&1
if [ $? -ne 0 ]; then
  echo "Build failed."
  exit 1
fi


cd ./fminer/ ; ./configure >/dev/null
if [ $? -ne 0 ]; then
  echo "Configuration failed."
  exit 1
fi
cd -
make -C fminer/ > /dev/null 2>&1
if [ $? -ne 0 ]; then
  echo "Build failed."
  exit 1
fi
