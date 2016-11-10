#!/bin/sh
#
# for use with the procedure in pandeia/doc/install/install.txt
#
# ddir = data directory
# sdir = source directory
# idir = install directory

cd $sdir

mkdir -p $idir/client
cp -r ui/client/. $idir/client/.

mkdir -p $idir/refdata

mkdir -p $ddir/results $ddir/log

[ -L $idir/results ] && rm -f $idir/results
[ -L $idir/log     ] && rm -f $idir/log 

mkdir -p $ddir/results
mkdir -p $ddir/log

[ -d $idir/results ] || ln -s $ddir/results $idir/results
[ -d $idir/log     ] || ln -s $ddir/log     $idir/log

mkdir -p $idir/db

cp $sdir/ui/ctrl* $idir/

cp $sdir/ui/ctrl_kick_it $idir/ctrl_kick_it

