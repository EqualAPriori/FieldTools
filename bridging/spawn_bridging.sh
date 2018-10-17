#!/bin/bash

bridgingdir=~/miktoarm/SCFT_Mikto/bridging/bridging-scripts
idir=`pwd`
#phiAmin=0.1
#phiAmax=0.8
calcminphase_flag=0
tau=0.5
chiN=30.000 #60.000
skip_if_bridging_exists=0
submit_flag=1
phiAmin=0.14 #0.14 #0.14
phiAmax=0.48 #0.64 #0.48
phiAdelta=0.02

#for i  in tau0.9_chiN30.000/phiA*; do
#for i  in tau0.5_chiN40.000/phiA*; do
#for i  in `seq 0.1 0.02 0.62`; do
for i  in `seq ${phiAmin} ${phiAdelta} ${phiAmax}`; do
    phiA=`echo $i | sed 's/0*$//g'`
    #wdir="tau${tau}_chiN${chiN}/phiA${phiA}"
    wdir="narms3_chiN${chiN}/phiA${phiA}"
    #phiA=`echo $wdir | cut -f 2 -d '/' | sed 's/[a-z,A-Z]//g' | xargs printf %f`
    #wdir="tau0.5_chiN60.000/phiA${phiA}"
    #wdir=$i
    if [ ! -e $wdir ]; then
        echo "Error! $wdir does not exist, cannot start a bridging calculation from here! Skipping..."
        continue
    fi
    echo Spawning bridging from $wdir 

    if [ $calcminphase_flag -eq 1 ]; then
        #$bridgingdir/get_minF_phase.py -d $wdir -o 
        ~/tools/phase-diagrams/extractF0.py -d $wdir --writemin
    fi

    minphase_file="$wdir/minphase.dat"
    if [ ! -e $minphase_file ]; then
        echo Skipping $wdir! $minphase_file does not exist! 
        continue
    fi
    
    minphase=`cat $minphase_file`
    if [ $minphase == "DIS" ]; then
       echo "Skipping $wdir. minphase == DIS"
       continue
    elif [ $minphase == "GYR" ]; then
       echo "Skipping $wdir. minphase == GYR"
       continue
    fi

    minphase_dir="${wdir}/${minphase}Phase"
    bridging_dir="$minphase_dir/bridging"
   

    if [ $skip_if_bridging_exists -eq 1 ] && [ -e $bridging_dir ]; then
        echo "$wdir/$minphase_dir/$bridging_dir already exists. Skipping!"
        continue
    else
        mkdir -p $bridging_dir
    fi

    cd $bridging_dir

    # deside whether to use GPU or CPU 
    if [ $minphase == "BCC" ]; then
        #submitfile="submitGPU.sh"
        submitfile="submit.sh"
    else
        submitfile="submit.sh"
    fi

    cmd="sed -e 's/__PHASE__/$minphase/g' 
             -e 's/__OUTFILE__/bridging.out/g'
             <$bridgingdir/$submitfile >submit.sh"
    eval $cmd
    
    if [ $submit_flag -eq 1 ]; then
        qsub submit.sh
    fi
    
    cd $idir
done
