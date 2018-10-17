#!/bin/bash

# ---------------------------------------------------------
# INPUT ARGS
# ---------------------------------------------------------
if [ $# -lt 2 ]; then
    echo "Usage: $0 <phase> <working path (with the SCFT outputs)> [density-threshold]"
    exit 1
fi

phase=$1 #"HEX"
mypath=$2
#runarg=$3
#region_to_analyze=${4}
density_threshold=${3:-0.50}
region_method="voronoi_burn"
#delta_density_threshold=0.05 # how much density threshold is automatically increased if too low
#max_density_threshold=0.5
#if [ $runarg == "all" ]; then
#    runa=1; runb=1; runc=1; rund=1;
#elif [ $runarg == "a" ]; then  runa=1; runb=0; runc=0; rund=0;
#elif [ $runarg == "b" ]; then  runa=0; runb=1; runc=0; rund=0;
#elif [ $runarg == "c" ]; then  runa=0; runb=0; runc=1; rund=0;
#elif [ $runarg == "d" ]; then  runa=0; runb=0; runc=0; rund=1; 
#elif [ $runarg == "cd" ]; then  runa=0; runb=0; runc=1; rund=1;
#else
#    echo "Error: Invalid runarg $runarg"
#    exit 1
#fi


idir=`pwd`
polyftsdir="/home/lequieu/PolyFTS/bin/Release"
#polyftsdir_bridging="/home/lequieu/PolyFTS-bridging/bin/Debug/"
polyftsdir_bridging="/home/lequieu/PolyFTS-bridging/bin/Release/"
#bridging_scripts="/home/lequieu/PolyFTS-bridging/bin/Release/doc/tutorial/41.Bridging_MiktoStar/scripts/"
bridging_scripts="/home/lequieu/miktoarm/SCFT_Mikto/bridging/bridging-scripts/"

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

# check for ${phase}.out and operators.dat
outfile="$idir/$mypath/${phase}.out"
if [ ! -s $outfile ]; then
    echo "ERROR! file $outfile does not exist!"; exit
fi
operatorfile="$idir/$mypath/operators.dat"
if [ ! -s $operatorfile ]; then
    echo "ERROR! file $operatorfile does not exist!"; exit
fi
fieldsfile="$idir/$mypath/fields.bin"
if [ ! -s $fieldsfile ]; then
    echo "ERROR! file $fieldsfile does not exist!"; exit
fi





# ---------------------------------------------
# Set some initial parameters
# ---------------------------------------------
if [ $phase == "SIGMA" ]; then ndomains=30
elif [ $phase == "A15" ]; then ndomains=8
elif [ $phase == "C15" ]; then ndomains=24
elif [ $phase == "C14" ]; then ndomains=12

elif [ $phase == "LAM" ] || [ $phase == "HEX" ] || [ $phase == "BCC" ]; then
    echo "Error! Invalid phase $phase. Are you sure you dont want to use the 'clasical' phases bridiging-wrapper.sh?"
    exit 1
else
    echo "Error! unknown phase $phase!"
    exit
fi

# ---------------------------------------------
# d.bridging
# ---------------------------------------------
ddir="d.bridging"
if [ 1 ]; then
    echo "Setting up $ddir..."
    mkdir -p $ddir
    cd $ddir
    cp $idir/$mypath/${phase}.in ${phase}.in.orig
    cp $idir/$mypath/fields.bin fields.in

    cmd="sed  
             -e '/NumBlocks/s/=.*/=\ 1/g' 
             -e '/NumTimeStepsPerBlock/s/=.*/=\ 1/g' 
             -e '/operators/a \ \ \ \ CalcBridging\ =\ true' 
             <${phase}.in.orig >${phase}.in"
    echo $cmd
    eval $cmd

     # cubic phases
    if [ $phase == "BCC" ] || [ $phase == "GYR" ] || [ $phase == "A15" ] || [ $phase == "FCC" ] || [ $phase == "C15" ]; then
        boxl=`tail -n 1 $operatorfile | awk '{print $10}'`
        sed  -i -e "/cellscaling/s/=.*/=\ $boxl/g" ${phase}.in
    #tetragonal phases
    elif [ $phase == "SIGMA" ]; then
        boxl=`tail -n 1 $operatorfile | awk '{print $10}'`
        zrel=`tail -n 1 $operatorfile | awk '{print $18/$10}'`
        sed  -i -e "/cellscaling/s/=.*/=\ $boxl/g" -e "/celllengths/s/=.*/=\ 1.0\ 1.0 $zrel/g" ${phase}.in
    # 3D hexagonal phases
    elif [ $phase == "C14" ]; then
        boxl=`tail -n 1 $operatorfile | awk '{print sqrt($10*$10+$11*$11)}'`
        zrel=`tail -n 1 $operatorfile | awk '{print $18/sqrt($10*$10+$11*$11)}'`
        sed  -i -e "/cellscaling/s/=.*/=\ $boxl/g" -e "/celllengths/s/=.*/=\ 1.0\ 1.0 $zrel/g" ${phase}.in
    fi


    domainsfile="domains_bridging.dat"
    rm $domainsfile
    for idomain in `seq 1 $ndomains`; do

        echo "Starting PolyFTS to compute bridging (idomain $idomain, density_threshold $density_threshold)"
        
        # create wdir and necessary input files
        workingdir=$idomain
        mkdir -p $workingdir
        cd $workingdir
        if [ -e fields.in ] ; then rm fields.in; fi
        ln -s ../fields.in .

        filename="bridging_params.in"
        echo "relative_density_threshold = $density_threshold" > $filename
        echo "region_to_analyze = $idomain" >> $filename
        echo "region_method = ${region_method}" >> $filename

        $polyftsdir_bridging/PolyFTS.x ../${phase}.in > ${phase}.out

        # if BridgingOperator.dat then it was a successful bridging
        if [ -s "BridgingOperator.dat" ]; then 
            # now integrate BridgingOperator.dat to get bridging fractions (and Pm)
            ${bridging_scripts}/integrate_Pmr.py > Pm.log
        fi 

        sumPtilda=`grep 'Ptilda' Pm.log | cut -f 2 -d ' '`
        sumPhat=`grep 'Phat' Pm.log | cut -f 2 -d ' '`
        echo $idomain $sumPtilda $sumPhat >> ../$domainsfile

        cd ..

    done
   
    cd $idir
fi

