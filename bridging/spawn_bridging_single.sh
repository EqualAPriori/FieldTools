#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Usage: $0 <working directory to start bridging simulation in> <is FK>"
    exit 1
fi
#wdir=tau0.9_chiN40.000/phiA0.2/BCCPhase
wdir=$1
isFKphase=$2
if [ ! -e $wdir ]; then
    echo "Error! $wdir does not exist!"
    exit 1
fi

bridgingscripts=~/miktoarm/SCFT_Mikto/bridging/bridging-scripts
if [ ! -e $bridgingscripts ]; then
    echo "Error! $bridgingscripts does not exist!"
    exit 1
fi

idir=`pwd`

bridging_dir="$wdir/bridging"
phase=`echo $wdir | sed 's:/:\n:g' | grep Phase | sed 's/Phase//g'`

mkdir -p $bridging_dir
cd $bridging_dir

# deside whether to use GPU or CPU 
    #submitfile="submitGPU.sh"
if [ $isFKphase -eq 1 ];then
    submitfile="submit-FK.sh"
else
    submitfile="submit.sh"
fi

cmd="sed -e 's/__PHASE__/${phase}/g' 
         -e 's/__OUTFILE__/bridging.out/g'
         <$bridgingscripts/$submitfile >submit.sh"
eval $cmd

qsub submit.sh

cd $idir
