#!/bin/bash
# Try running this file for a complete build

make -C liblast/ clean >/dev/null 2>&1
make -C libbbrc/ clean >/dev/null 2>&1
make -C fminer/ clean >/dev/null 2>&1

echo "Configuring LibLAST..."
cd ./liblast/ ; ./configure >> build.log 2>&1
if [ $? -ne 0 ]; then
  echo "Configuration failed. Please examine 'liblast/build.log'."
  exit 1
fi
cd -
echo "Building LibLAST..."
if ! make -C liblast/ >> liblast/build.log 2>&1; then
  echo "Build failed. Please examine 'liblast/build.log'."
  exit 1
fi
echo "Building LibLAST for Ruby..."
if ! make -C liblast/ ruby >> liblast/build.log 2>&1; then
  echo "Build failed. Please examine 'liblast/build.log'."
  exit 1
fi

echo "Configuring LibBBRC..."
cd ./libbbrc/ ; ./configure >> build.log 2>&1
if [ $? -ne 0 ]; then
  echo "Configuration failed. Please examine 'libbbrc/build.log'."
  exit 1
fi
cd -
echo "Building LibBBRC..."
make -C libbbrc/ >> libbbrc/build.log 2>&1
if [ $? -ne 0 ]; then
  echo "Build failed. Please examine 'libbbrc/build.log'."
  exit 1
fi
echo "Building LibBBRC for Ruby..."
make -C libbbrc/ ruby >> libbbrc/build.log 2>&1
if [ $? -ne 0 ]; then
  echo "Build failed. Please examine 'libbbrc/build.log'."
  exit 1
fi


echo "Configuring Fminer..."
cd ./fminer/ ; ./configure >> build.log 2>&1
if [ $? -ne 0 ]; then
  echo "Configuration failed. Please examine 'fminer/build.log'."
  exit 1
fi
cd -
echo "Building Fminer..."
make -C fminer/ >> fminer/build.log 2>&1
if [ $? -ne 0 ]; then
  echo "Build failed. Please examine 'build.log'."
  exit 1
fi
