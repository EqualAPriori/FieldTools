#!/bin/bash
idir=`pwd`

# ---------------------------------------------------------
# INPUT ARGS
# ---------------------------------------------------------
if [ $# -lt 2 ]; then
    echo "Usage: $0 <working path (with the SCFT outputs)> <runarg [all|a|b|c|d]> [region-to-analyze] [density-threshold]"
    exit 1
fi

mypath=$1
runarg=$2
region_to_analyze=${3:-2}
density_threshold=${4:-0.1}
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
polyftsdir="/home/lequieu/PolyFTS/bin/Release/"
polyftsdir_bridging="/home/lequieu/PolyFTS-bridging/bin/Debug/"
phase="LAM"
nx=3
ny=1
nz=1
dim=1
# ---------------------------------------------------------


if [ $dim -eq 1 ]; then ny=1; nz=1; fi
if [ $dim -eq 2 ]; then nz=1; fi

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
             -e '/ReadInputFields/s/=.*/=\ Hfields\nInputFieldsFile = fields.in/g' 
             -e '/OutputFields/a \ \ \ \ OutputFormattedFields\ =\ true' 
             <${phase}.in.orig >${phase}.in"
    echo $cmd
    eval $cmd

    $polyftsdir/PolyFTS.x ${phase}.in 

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
    ~/PolyFTS-bridging/tools/PeriodicallyReplicateField.pl fields.dat $nx $ny $nz
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
    lxnew=`echo $lx*$nx | bc -l`
    npwx=`grep npw ${phase}.in.orig | sed 's/.*=//g' | awk '{print $1}'`
    npwxnew=`echo $npwx*$nx | bc -l`
    if [ $dim -ge 2 ]; then
        ly=`grep celllengths ${phase}.in.orig | sed 's/.*=//g' | awk '{print $2}'`
        lynew=`echo $ly*$ny | bc -l`
        npwy=`grep npw ${phase}.in.orig | sed 's/.*=//g' | awk '{print $2}'`
        npwynew=`echo $npwy*$ny | bc -l`
    fi
    if [ $dim -eq 3 ]; then
        lz=`grep celllengths ${phase}.in.orig | sed 's/.*=//g' | awk '{print $3}'`
        lznew=`echo $lz*$nz | bc -l`
        npwz=`grep npw ${phase}.in.orig | sed 's/.*=//g' | awk '{print $3}'`
        npwznew=`echo $npwz*$nz | bc -l`
    fi
    cmd="sed  -e '/celllengths/s/=.*/=\ $lxnew/g'
              -e '/npw/s/=.*/=\ $npwxnew/g'
              -e '/NumBlocks/s/=.*/=\ 1000/g' 
              -e '/NumTimeStepsPerBlock/s/=.*/=\ 100/g' 
              -e '/FieldOutputSpace/s/=.*/=\ both/g'
              < ${phase}.in.orig > ${phase}.in"
    echo $cmd
    eval $cmd

    $polyftsdir/PolyFTS.x ${phase}.in

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
    
    filename="bridging_params.in"
    if [ -e $filename ]; then rm $filename; fi 

    echo "relative_density_threshold = $density_threshold" >> $filename
    echo "region_to_analyze = $region_to_analyze" >> $filename

    $polyftsdir_bridging/PolyFTS.x ${phase}.in

    cd $idir
fi
