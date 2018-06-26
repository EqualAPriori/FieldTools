#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 <status to resubmit [0,1,2,3] >"
    echo "   0 - still running (typically hit max wall time or seg faulted) "
    echo "   1 - divergent trajectory "
    echo "   2 - SCFT converged "
    echo "   3 - ran to max timesteps without converging "
    exit
fi


idir=`pwd`
factor=10
status2resub=$1

function change_cell_lengths() {
    if [ $# -ne 2 ]; then
        echo "Usage: $0 <phase> <working directory>"
        exit
    fi
    phase=$1
    wdir=$2
    
    op_filename="$wdir/operators.dat"
    if [ ! -e $op_filename ]; then
        pwd
        echo "Error! $op_filename does not exist!"; exit
    fi

    #backup previous PHASE.in file
    in_filename=$wdir/${phase}.in
    index=`find $wdir -maxdepth 1 -name "${phase}.in.orig*" -type f | wc -l`
    in0_filename=$wdir/${phase}.in.orig.${index}
    if [ ! -e $in_filename ]; then
        echo "Error! $in_filename does not exist!"; exit
    fi
    echo $in0_filename
    if [ -e ${in0_filename} ]; then
        echo "Error! ${in0_filename} already exists, looks like this WDIR has already been resubmitted! (or was partially submitted!)";
        exit
    fi
    cp ${in_filename} ${in0_filename}
    
    # cubic phases
    if [ $phase == "BCC" ] || [ $phase == "GYR" ] || [ $phase == "A15" ] || [ $phase == "FCC" ] || [ $phase == "C15" ]; then
        boxl=`tail -n 1 $op_filename | awk '{print $10}'`
        sed  -i -e "/cellscaling/s/=.*/=\ $boxl/g" ${in_filename}
    #tetragonal phases
    elif [ $phase == "SIGMA" ]; then
        boxl=`tail -n 1 $op_filename | awk '{print $10}'`
        zrel=`tail -n 1 $op_filename | awk '{print $18/$10}'`
        sed  -i -e "/cellscaling/s/=.*/=\ $boxl/g" -e "/celllengths/s/=.*/=\ 1.0\ 1.0 $zrel/g" ${in_filename}
    # 2D hexagonal phases
    elif [ $phase == "HEX" ]; then
        boxl=`tail -n 1 $op_filename | awk '{print sqrt($7*$7+$8*$8)}'`
        sed  -i -e "/cellscaling/s/=.*/=\ $boxl/g" ${in_filename}
    # 3D hexagonal phases
    elif [ $phase == "C14" ]; then
        boxl=`tail -n 1 $op_filename | awk '{print sqrt($10*$10+$11*$11)}'`
        zrel=`tail -n 1 $op_filename | awk '{print $18/sqrt($10*$10+$11*$11)}'`
        sed  -i -e "/cellscaling/s/=.*/=\ $boxl/g" -e "/celllengths/s/=.*/=\ 1.0\ 1.0 $zrel/g" ${in_filename}
    # 1D phases
    elif [ $phase == "LAM" ] || [ $phase == "DIS" ] ; then
        boxl=`tail -n 1 $op_filename | awk '{print $5}'`
        sed  -i -e "/celllengths/s/=.*/=\ $boxl/g" ${in_filename}
    # unknown phase
    else
        echo "Error! \"$phase\" is an unknown phase!"
        exit
    fi
}

# ---------------------------------------------------------------------
#                        MAIN METHOD
# ---------------------------------------------------------------------
for i in ./tau0.9_chiN*/phiA*/*Phase/ ./chiN*/phiA*/*Phase/; do
    if [ ! -s $i/STATUS ]; then 
        echo $i/STATUS does not exist. Skipping!
        continue
    fi

    status=`cat $i/STATUS`
    if [[ $status -eq $status2resub ]]; then

        cd $i
        
        # get phase name
        phase=`echo "$i" | awk 'BEGIN {FS="/"}{print $(NF-1)}' | sed 's/Phase//g'`
    
        # check runtime and compare to walltime
        #runtimeline=`grep 'TOTAL Runtime' ${phase}.out | tail -n1`
        #runtime=`echo $runtimeline | awk '{print $(NF-1)}'` #runtime in seconds
        #runtime=`echo $runtime/3600. | bc -l` # runtime in hours

        #walltime=`grep 'walltime' submit.sh | head -n1 | awk 'BEGIN{FS="="}{print $NF}' | cut -f 1 -d ':'`
        #newruntime=`echo $runtime*$factor | bc -l`
        
        # if newruntime is bigger than walltime, then update walltime
        #status=`echo $newruntime'>'$walltime | bc -l`
        #if [ $status ]; then
        #    newwalltime=`echo $newruntime*2/1 | bc `
        #    if [ $newwalltime -eq 0 ]; then newwalltime=1; fi
        #    newwalltime=24
        #    echo "     Changing walltime to $newwalltime:00:00"
        #    cmd="sed -i '/PBS.*walltime/s/=.*/=$newwalltime:00:00/g' submit.sh"
        #    eval $cmd
        #fi

        if [[ $status -eq 0 ]] || [[ $status -eq 3 ]]; then
            echo "Resubmitting $i"
            
            # change fields.in
            index=`find . -maxdepth 1 -name 'fields.in.orig*'  -type f | wc -l`
            cp fields.in fields.in.orig.$index
            cp fields.bin fields.in

            #change cell lengths
            change_cell_lengths $phase .
            
            # change walltime
            newwalltime=24
            echo "     Changing walltime to $newwalltime:00:00"
            cmd="sed -i '/PBS.*walltime/s/=.*/=$newwalltime:00:00/g' submit.sh"
            eval $cmd


             change NumBlocks
            numblocks=`grep 'NumBlocks' ${phase}.in | awk '{FS="="}{print $NF}'`
            newnumblocks=`echo $numblocks*$factor | bc`
            newnumblocks=10000
            cmd="sed -i '/NumBlocks/s/=.*/=\ $newnumblocks/g' ${phase}.in"
            echo $cmd
            eval $cmd
    
            qsub submit.sh 


        elif [[ $status -eq 1 ]]; then

            newnumblocks=1000
            newDT=0.1
            cmd="sed -i -e 's///g' 
                        -e '/NumBlocks/s/=.*/=\ $newnumblocks/g'
                        -e '/TimeStepDT/s/=.*/=\ $newDT/g'
                            ${phase}.in"
            echo $cmd
            eval $cmd

            newwalltime=24
            cmd="sed -i '/PBS.*walltime/s/=.*/=$newwalltime:00:00/g' submit.sh"
            echo $cmd
            eval $cmd
        
            qsub submit.sh 


        fi 

        cd $idir
    fi
done


