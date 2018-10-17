#!/bin/bash

# ---------------------------------------------------------
# INPUT ARGS
# ---------------------------------------------------------
if [ $# -lt 2 ]; then
    echo "Usage: $0 <phase> <working path (with the SCFT outputs)> <runarg [all|a|b|c|d]> [region-to-analyze] [density-threshold]"
    exit 1
fi

phase=$1 #"HEX"
mypath=$2
runarg=$3
region_to_analyze=${4}
density_threshold=${5:-0.50}
delta_density_threshold=0.05 # how much density threshold is automatically increased if too low
max_density_threshold=0.50
region_method="voronoi_center"
if [ $runarg == "all" ]; then
    runa=1; runb=1; runc=1; rund=1;
elif [ $runarg == "a" ]; then  runa=1; runb=0; runc=0; rund=0;
elif [ $runarg == "b" ]; then  runa=0; runb=1; runc=0; rund=0;
elif [ $runarg == "c" ]; then  runa=0; runb=0; runc=1; rund=0;
elif [ $runarg == "d" ]; then  runa=0; runb=0; runc=0; rund=1; 
elif [ $runarg == "cd" ]; then  runa=0; runb=0; runc=1; rund=1;
else
    echo "Error: Invalid runarg $runarg"
    exit 1
fi


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



# ---------------------------------------------
# Set some initial parameters
# ---------------------------------------------
if [ $phase == "LAM" ]; then
    dim=1
    nx=3; ny=1; nz=1
    if [ -z $region_to_analyze ]; then
        region_to_analyze=2
    fi
elif [ $phase == "HEX" ]; then
    dim=2
    nx=3; ny=2; nz=1
    if [ -z $region_to_analyze ]; then
        region_to_analyze=9
    fi
elif [ $phase == "BCC" ]; then
    dim=3
    nx=2; ny=2; nz=2
    if [ -z $region_to_analyze ]; then
        region_to_analyze=12
    fi
else
    echo "Error! unknown phase $phase!"
    exit
fi
echo region_to_analyze: $region_to_analyze density_threshold: $density_threshold

# get boxtype
if [ $phase == "LAM" ] || [ $phase == "BCC" ] ; then
    orthorhombicbox=1
elif [ $phase == "HEX" ]; then
    if [ $dim == 2 ]; then
        offdiagzero=`tail -n 1 $operatorfile | awk '{print sqrt($8*$8)<1e-3}'` 
        echo offdiag $offdiagzero
        if [ $offdiagzero == 1 ]; then
            orthorhombicbox=1
        else
            orthorhombicbox=0
        fi
    fi
fi
echo runa: $runa orthorhombicbox: $orthorhombicbox
# ---------------------------------------------
#  a (depends on if box was orthorhombic or not)
# ---------------------------------------------

if [ $runa == 1 ] && [ $orthorhombicbox == 1 ]; then
    # ---------------------------------------------
    # a.format-fields
    # ---------------------------------------------
    adir="a.format-fields"
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
    boxl=`tail -n 1 $operatorfile | awk '{print $10}'`

    # change boxl in ${phase}.in, and change to orthorhombic
    cmd="sed -e 's///g'  
             -e '/cellscaling/s/=.*/=\ ${boxl}/g'       
             -e '/NumBlocks/s/=.*/=\ 1/g' 
             -e '/NumTimeStepsPerBlock/s/=.*/=\ 1/g' 
             -e '/OutputFields/a \ \ \ \ OutputFormattedFields\ =\ true' 
             <${phase}.in.orig >${phase}.in"
    echo $cmd
    eval $cmd

    #LAM doesn't read fields so need to change ${phase}.in
    if [ $phase == "LAM" ]; then 
        cmd="sed -i '/ReadInputFields/s/=.*/=\ Hfields\nInputFieldsFile = fields.in/g' ${phase}.in"
        echo $cmd
        eval $cmd
    fi

    $polyftsdir/PolyFTS.x ${phase}.in 

    cd $idir

elif [ $runa -eq 1 ]; then
    # ---------------------------------------------
    # a.make-orthorhombic-seed
    # ---------------------------------------------
    #  need to make orthorhombic box in order to replicate fields in "step b"

    adir="a.make-orthorhombic-seed"
    echo "Setting up $adir..."
    mkdir -p $adir
    cd $adir

    infileorig=$idir/$mypath/${phase}.in 
    if [ ! -e $infileorig ]; then
        echo "Error! file $infileorig does not exist!"
        exit 1
    fi
    
    cp -r $infileorig ${phase}.in.orig

    #get final box length (note hexagonal)
    boxl=`tail -n 1 $operatorfile | awk '{print sqrt($7*$7+$8*$8)}'` #get absolute length since hexagonal box

    # change boxl in ${phase}.in, and change to orthorhombic
    cmd="sed -e 's///g'  
             -e '/cellscaling/s/=.*/=\ ${boxl}/g'       
             -e '/celllengths/s/=.*/=\ 1.0\ 1.73/g' 
             -e '/npw/s/=.*/=\ 32\ 48/g' 
             -e '/cellangles/s/=.*/=\ 90.0/g' 
             -e '/TimeStepDT/s/=.*/=\ 1.00/g' 
             -e '/spacegroupname/d ' 
             -e '/symmetrize/s/=.*/=\ off/g' 
             -e '/NumBlocks/s/=.*/=\ 500/g' 
             -e '/OutputFields/a \ \ \ \ OutputFormattedFields\ =\ true' 
             -e '/ReadInputFields/s/=.*/=\ no\n\tInitField1{\n\tInitType=URNG\n\t}\n\tInitField2{\n\tInitType=URNG\n\t}\n/g' 
             <${phase}.in.orig >${phase}.in"
    echo $cmd
    eval $cmd
    
    #loop PolyFTS call until SCFT converges
    done=0; count=1
    while [ $done -eq 0 ]; do
        echo "Starting PolyFTS to make orthorhombic seed (attempt $count)"
        $polyftsdir/PolyFTS.x ${phase}.in > ${phase}.out

        status=`cat STATUS`
        if [ $status -eq 2 ]; then done=1; fi # SCFT converged
        
        #of not done, copy attempt files into directory
        if [ $done -eq 0 ]; then
            mkdir -p attempt$count
            cp * attempt$count
        fi

        if [ $count -gt 20 ]; then
            echo "Error! Too many attempts making orthorhombic seed. Something looks wrong!"
            exit
        fi



        let count+=1
    done
    
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

    #change npw and cell lengths
    if [ $dim -eq 1 ]; then 
        cmd="sed  -e '/celllengths/s/=.*/=\ $lxnew/g'
                  -e '/npw/s/=.*/=\ $npwxnew/g'
                  < ${phase}.in.orig > ${phase}.in"
    elif [ $dim -eq 2 ]; then 
        cmd="sed  -e '/celllengths/s/=.*/=\ $lxnew\ $lynew/g'
                  -e '/npw/s/=.*/=\ $npwxnew\ $npwynew/g'
                  < ${phase}.in.orig > ${phase}.in"

    elif [ $dim -eq 3 ]; then 
        cmd="sed  -e '/celllengths/s/=.*/=\ $lxnew\ $lynew\ $lznew/g'
                  -e '/npw/s/=.*/=\ $npwxnew\ $npwynew\ $npwznew/g'
                  < ${phase}.in.orig > ${phase}.in"
    fi
    echo $cmd
    eval $cmd
    
    
    # now change runtime parameters, using new input file!
    cmd="sed  -i -e '/NumBlocks/s/=.*/=\ 500/g' 
                 -e '/NumTimeStepsPerBlock/s/=.*/=\ 100/g' 
                 -e '/FieldOutputSpace/s/=.*/=\ both/g'
               ${phase}.in"
    echo $cmd
    eval $cmd

    #if not orthorhombic box, need to turn on field reading again
    if [ $orthorhombicbox == 0 ]; then
        cmd="sed  -i '/ReadInputFields/s/=.*/=\ yes/g' ${phase}.in"
        echo $cmd 
        eval $cmd 

    fi

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
        if [ $count -gt 20 ]; then
            echo "Error! Too many attempts at large cell SCFT. Something looks wrong!"
            exit
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
    #cp $idir/$cdir/fields.dat fields.in
    cp $idir/$cdir/fields_k.dat fields.in

    cmd="sed -e 's///g'  
             -e '/NumBlocks/s/=.*/=\ 1/g' 
             -e '/NumTimeStepsPerBlock/s/=.*/=\ 1/g' 
             -e '/operators/a \ \ \ \ CalcBridging\ =\ true' 
             <${phase}.in.orig >${phase}.in"
    echo $cmd
    eval $cmd
    
    if [ $phase == "LAM" ]; then
        cmd="sed -i '/npw/s/=.*/= 256/g' ${phase}.in"
        echo $cmd
        eval $cmd
    fi
    

    count=0
    #threshfile="threshold_vs_bfraction.dat"
    #rm $threshfile

    for threshold in `seq -w $density_threshold $delta_density_threshold $max_density_threshold`; do

        echo "Starting PolyFTS to compute bridging (attempt $count, density_threshold $threshold)"
        
        # create wdir and necessary input files
        #workingdir="voronoi-$threshold"
        workingdir="voronoi"
        mkdir -p $workingdir
        cd $workingdir
        if [ -e fields.in ] ; then rm fields.in; fi
        ln -s ../fields.in .

        filename="bridging_params.in"
        echo "relative_density_threshold = $threshold" > $filename
        echo "region_to_analyze = $region_to_analyze" >> $filename
        echo "region_method = ${region_method}" >> $filename

        $polyftsdir_bridging/PolyFTS.x ../${phase}.in > ${phase}.out

        # if BridgingOperator.dat then it was a successful bridging
        if [ -s "BridgingOperator.dat" ]; then 
            # now integrate BridgingOperator.dat to get bridging fractions (and Pm)
            ${bridging_scripts}/integrate_Pmr.py > Pm.log

            let count+=1
            if [ $count -eq 1 ];then
                fitstart=$threshold
            elif [ $count -eq 5 ];then
                fitend=$threshold
            fi
        fi 

        #sumPtilda=`grep 'sumPtilda' Pm.log | cut -f 2 -d ' '`
        #sumPhat=`grep 'sumPhat' Pm.log | cut -f 2 -d ' '`
        #echo $threshold $sumPtilda $sumPhat >> ../$threshfile

        cd ..

    done

    
#    # now we should have a series of BridgingFractions for different thresholds
#    # now need to fit a line and extract zero-threshold intercept
#    if [ $count -lt 5 ]; then
#        echo "Error! Less that 5 valid points for fit"
#        exit 1    
#    fi
#    echo fitstart: $fitstart, fitend: $fitend
#
##really shouldn't use gnuplot for the fit, but its a quick and dirty solution
#gnuplot << EOF
#f(x)=a*x+b
#fit [$fitstart:$fitend] f(x) './$threshfile' u 1:2 via a,b
#set print "bridging_fraction.dat"
#print b
#EOF
#    
#    echo "Bridging Fraction from fitting line is `cat bridging_fraction.dat`"

    cd $idir
fi

