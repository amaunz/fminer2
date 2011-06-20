#!/bin/sh
# vim:et:ft=sh:sts=2:sw=2
#
# Copyright 2008 Kate Ward. All Rights Reserved.
# Released under the LGPL (GNU Lesser General Public License)
#
# Author: kate.ward@forestent.com (Kate Ward)
#
# Example unit test for the mkdir command.
#
# There are times when an existing shell script needs to be tested. In this
# example, we will test several aspects of the the mkdir command, but the
# techniques could be used for any existing shell script.

#-----------------------------------------------------------------------------
# suite tests
#

testBbrcClassification()
{
  mkdir $testdir >/dev/null 2>&1
  $fminer $libbrc $hamster>$testdir/tmp1 2>$testdir/tmp1e
  h=`md5sum $testdir/tmp1 | sed 's/\s.*//g'`
  he=`md5sum $testdir/tmp1e | sed 's/\s.*//g'`
  assertEquals "testBbrcClassification" "$h" "3b8c0f319d38be021add34b62fda2b44"
}

testBbrcFsmClassification()
{
  mkdir $testdir >/dev/null 2>&1
  $fminer $libbrc $fsmargs $hamster>$testdir/tmp2 2>$testdir/tmp2e
  h=`md5sum $testdir/tmp2 | sed 's/\s.*//g'`
  he=`md5sum $testdir/tmp2e | sed 's/\s.*//g'`
  assertEquals "testBbrcClassification" "$h" "576522f98236969193b8f661309ccefc"
}

testBbrcFsmP0Classification()
{
  mkdir $testdir >/dev/null 2>&1
  $fminer $libbrc $fsmargs $p0args $hamster>$testdir/tmp3 2>$testdir/tmp3e
  h=`md5sum $testdir/tmp3 | sed 's/\s.*//g'`
  he=`md5sum $testdir/tmp3e | sed 's/\s.*//g'`
  assertEquals "testBbrcClassification" "$h" "cf5a20167243312554c26127af22954e"
}

#-----------------------------------------------------------------------------
# suite functions
#

oneTimeSetUp()
{
  fminer='./fminer'  # save command name in variable to make future changes easy
  libbrc='../libbbrc/libbbrc.so'
  liblast='../liblast/liblast.so'
  fsmargs='-d -b'
  p0args='-p0'
  hamster='../libbbrc/test/hamster_carcinogenicity.smi ../libbbrc/test/hamster_carcinogenicity.class'
  epafhm='../libbbrc/test/EPAFHM.smi ../libbbrc/test/EPAFHM.act'

  md5='md5sum | awk -F " " "{print $1}"'
  testdir="./test/tmp"
  unset FMINER_LAZAR
  export FMINER_SMARTS
}

tearDown()
{
  rm -fr "${testdir}"
}

# load and run shUnit2
[ -n "${ZSH_VERSION:-}" ] && SHUNIT_PARENT=$0
. "test/shunit/src/shunit2"
