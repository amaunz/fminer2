#!/bin/sh
# suite tests
#

testBbrcClassification()
{
  mkdir $testdir >/dev/null 2>&1
  $fminer $libbrc $hamster>$testdir/tmp1 2>$testdir/tmp1e
  h=`md5sum $testdir/tmp1 | sed 's/\s.*//g'`
  he=`md5sum $testdir/tmp1e | sed 's/\s.*//g'`
  if [ "`uname -m`" = "x86_64" ]; then
    assertEquals "testBbrcClassification" "$h" "3b8c0f319d38be021add34b62fda2b44"
  else
    assertEquals "testBbrcClassification" "$h" "599fb20d22d88226377047c15bfb9ca8"
  fi
}

testBbrcFsmClassification()
{
  mkdir $testdir >/dev/null 2>&1
  $fminer $libbrc $fsmargs $hamster>$testdir/tmp2 2>$testdir/tmp2e
  h=`md5sum $testdir/tmp2 | sed 's/\s.*//g'`
  he=`md5sum $testdir/tmp2e | sed 's/\s.*//g'`
  if [ "`uname -m`" = "x86_64" ]; then
    assertEquals "testBbrcClassification" "$h" "576522f98236969193b8f661309ccefc"
  else
    assertEquals "testBbrcClassification" "$h" "1ce56e91303503e6fd9cb5ecf1dc14fb"
  fi
}

testBbrcFsmP0Classification()
{
  mkdir $testdir >/dev/null 2>&1
  $fminer $libbrc $fsmargs $p0args $hamster>$testdir/tmp3 2>$testdir/tmp3e
  h=`md5sum $testdir/tmp3 | sed 's/\s.*//g'`
  he=`md5sum $testdir/tmp3e | sed 's/\s.*//g'`
  if [ "`uname -m`" = "x86_64" ]; then
    assertEquals "testBbrcClassification" "$h" "cf5a20167243312554c26127af22954e"
  else
    assertEquals "testBbrcClassification" "$h" "aad799b68060d8a6262230aae954f44c"
  fi
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
  echo -n ""
  #rm -fr "${testdir}"
}

# load and run shUnit2
[ -n "${ZSH_VERSION:-}" ] && SHUNIT_PARENT=$0
. "test/shunit/src/shunit2"
