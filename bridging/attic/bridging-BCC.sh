#!/bin/bash

# ---------------------------------------------------------
# INPUT ARGS
# ---------------------------------------------------------
if [ $# -lt 2 ]; then
    echo "Usage: $0 <working path (with the SCFT outputs)> <runarg [all|a|b|c|d]> [region-to-analyze] [density-threshold]"
    exit 1
fi

mypath=$1
runarg=$2
region_to_analyze=${3:-12}
density_threshold=${4:-0.05}
delta_density_threshold=0.05 # how much density threshold is automatically increased if too low
if [ $runarg == "all" ]; then
    runa=1; runb=1; runc=1; rund=1;
elif [ $runarg == "a" ]; then  runa=1; runb=0; runc=0; rund=0;
elif [ $runarg == "b" ]; then  runa=0; runb=1; runc=0; rund=0;
elif [ $runarg == "c" ]; then  runa=0; runb=0; runc=1; rund=0;
elif [ $runarg == "d" ]; then  runa=0; runb=0; runc=0; rund=1;
else
    echo "Error: Invalid runarg $runarg"
    exit 1
fi

echo region_to_analyze: $region_to_analyze density_threshold: $density_threshold

idir=`pwd`
polyftsdir="/home/lequieu/PolyFTS/bin/Release"
polyftsdir_bridging="/home/lequieu/PolyFTS-bridging/bin/Debug/"
phase="BCC"
nx=2
ny=2
nz=2
# ---------------------------------------------------------


# ---------------------------------------------
# Initial Checks
# ---------------------------------------------

# check that previous simulation converged
statusfile="$idir/$mypath/STATUS"
if [ ! -s $statusfile ];then
    echo "ERROR! file $statusfile does not exist!"
    exit
fi
status=`cat $statusfile`
if [ $status -ne 2 ]; then
    echo "ERROR! SCFT run in $idir/$mypath is not converged! Cannot use to calculate bridging!"
    exit
fi

# ---------------------------------------------
# a.format-fields
# ---------------------------------------------
adir="a.format-fields"
if [ $runa -eq 1 ]; then
    echo "Setting up $adir..."
    mkdir -p $adir
    cd $adir

    infileorig=$idir/$mypath/${phase}.in 
    if [ ! -e $infileorig ]; then
        echo "Error! file $infileorig does not exist!"
        exit 1
    fi
    
    cp -r $infileorig ${phase}.in.orig
    cp $idir/$mypath/fields.bin fields.in
    
    #get final box length
    boxl=`tail -n 1 $idir/$mypath/operators.dat | awk '{print $10}'`
    
    # change boxl in ${phase}.in, and change to orthorhombic
    cmd="sed -e 's///g'  
             -e '/cellscaling/s/=.*/=\ ${boxl}/g'       
             -e '/NumBlocks/s/=.*/=\ 1/g' 
             -e '/NumTimeStepsPerBlock/s/=.*/=\ 1/g' 
             -e '/OutputFields/a \ \ \ \ OutputFormattedFields\ =\ true' 
             <${phase}.in.orig >${phase}.in"
    echo $cmd
    eval $cmd
    
    $polyftsdir/PolyFTS.x ${phase}.in > ${phase}.out
    
    cd $idir
fi

# ---------------------------------------------
# b.replicate-fields
# ---------------------------------------------
bdir="b.replicate-fields"
if [ $runb -eq 1 ]; then
    echo "Setting up $bdir..."
    mkdir -p $bdir
    cd $bdir
    cp $idir/$adir/fields.dat .
    #~/PolyFTS-bridging/tools/PeriodicallyReplicateField.pl fields.dat $nx $ny $nz
    ${polyftsdir_bridging}/../../tools/PeriodicallyReplicateField.pl fields.dat $nx $ny $nz
    cd $idir
fi

# ---------------------------------------------
# c.scft-replicate-fields
# ---------------------------------------------
cdir="c.scft-replicated-fields"
if [ $runc -eq 1 ]; then
    echo "Setting up $cdir..."
    mkdir -p $cdir
    cd $cdir
    cp $idir/$adir/${phase}.in ${phase}.in.orig
    cp $idir/$bdir/fields.dat_replicated.dat fields.in

    #change celllengths, npw, to correspond to the replicated fields
    lx=`grep celllengths ${phase}.in.orig | sed 's/.*=//g' | awk '{print $1}'`
    ly=`grep celllengths ${phase}.in.orig | sed 's/.*=//g' | awk '{print $2}'`
    lz=`grep celllengths ${phase}.in.orig | sed 's/.*=//g' | awk '{print $3}'`
    lxnew=`echo $lx*$nx | bc -l`
    lynew=`echo $ly*$ny | bc -l`
    lznew=`echo $lz*$nz | bc -l`
    npwx=`grep npw ${phase}.in.orig | sed 's/.*=//g' | awk '{print $1}'`
    npwy=`grep npw ${phase}.in.orig | sed 's/.*=//g' | awk '{print $2}'`
    npwz=`grep npw ${phase}.in.orig | sed 's/.*=//g' | awk '{print $3}'`
    npwxnew=`echo $npwx*$nx | bc -l`
    npwynew=`echo $npwy*$ny | bc -l`
    npwznew=`echo $npwz*$nz | bc -l`
    cmd="sed  -e '/celllengths/s/=.*/=\ $lxnew\ $lynew\ $lznew/g'
              -e '/npw/s/=.*/=\ $npwxnew\ $npwynew\ $npwznew/g'
              -e '/NumBlocks/s/=.*/=\ 100/g' 
              -e '/NumTimeStepsPerBlock/s/=.*/=\ 100/g' 
              -e '/FieldOutputSpace/s/=.*/=\ both/g'
              < ${phase}.in.orig > ${phase}.in"
    echo $cmd
    eval $cmd

    #$polyftsdir/PolyFTS.x ${phase}.in > ${phase}.out
    #loop PolyFTS call until SCFT converges
    done=0; count=1
    while [ $done -eq 0 ]; do
        echo "Starting PolyFTS to converge replicated fields (attempt $count)"
        $polyftsdir/PolyFTS.x ${phase}.in > ${phase}.out

        status=`cat STATUS`
        if [ $status -eq 2 ]; then done=1; fi # SCFT converged

        #of not done, copy attempt files into directory
        if [ $done -eq 0 ]; then
            mkdir -p attempt$count
            cp * attempt$count
        fi

        let count+=1
    done


    cd $idir
fi

# ---------------------------------------------
# d.bridging
# ---------------------------------------------
ddir="d.bridging"
if [ $rund -eq 1 ]; then
    echo "Setting up $ddir..."
    mkdir -p $ddir
    cd $ddir
    cp $idir/$cdir/${phase}.in ${phase}.in.orig
    cp $idir/$cdir/fields.dat fields.in

    cmd="sed -e 's///g'  
             -e '/NumBlocks/s/=.*/=\ 1/g' 
             -e '/NumTimeStepsPerBlock/s/=.*/=\ 1/g' 
             -e '/operators/a \ \ \ \ CalcBridging\ =\ true' 
             <${phase}.in.orig >${phase}.in"
    echo $cmd
    eval $cmd
    
    #filename="bridging_params.in"
    #rm $filename
    #echo "relative_density_threshold = $density_threshold" >> $filename
    #echo "region_to_analyze = $region_to_analyze" >> $filename

    #$polyftsdir_bridging/PolyFTS.x ${phase}.in

    done=0; count=1
    while [ $done -eq 0 ]; do

        echo "Starting PolyFTS to compute bridging (attempt $count, density_threshold $density_threshold)"

        filename="bridging_params.in"
        rm $filename
        echo "relative_density_threshold = $density_threshold" >> $filename
        echo "region_to_analyze = $region_to_analyze" >> $filename

        #$polyftsdir/PolyFTS.x ${phase}.in > ${phase}.out.${count}
        $polyftsdir_bridging/PolyFTS.x ${phase}.in > ${phase}.out

        # if BridgingOperator.dat then we're done
        if [ -s "BridgingOperator.dat" ]; then done=1; fi 

        #of not done, copy attempt files into directory
        if [ $done -eq 0 ]; then
            mkdir -p attempt$count
            cp * attempt$count
        fi

        density_threshold=`echo "$density_threshold + $delta_density_threshold" | bc -l`
        let count+=1
    done



    cd $idir
fi
