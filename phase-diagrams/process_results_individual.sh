#!/bin/bash

#Extract free-energy vs fA for all different phases and plot the results
set -o pipefail # Allow error codes to pass through pipes
IDIR=`pwd`

filename=F0_phases.dat
for dir1 in chiN* tau* eps*; do
  cd $dir1

  for dir2 in phiA*; do
    cd $dir2

    rm -f ${filename}

    # Write free-energy data
    for dir3 in *Phase; do
      phase=$dir3
      # ignore certain phases?
      #if [[ $dir3 == "C14Phase" || $dir3 == "C15Phase" ]]; then continue; fi
      if [ -e $dir3/STATUS ]; then
          status=`cat $dir3/STATUS`
          if [ $status == "1" ]; then
               echo "$dir1/$dir2/$dir3 is divergent....skipping!" 
               continue
          #elif [ $status == "3" ]; then
          #     echo "$dir1/$dir2/$dir3 SCFT is not converged....skipping!" 
          #     continue
          fi
      else
        echo "$dir1/$dir2/$dir3/STATUS does not exist...skipping!"
        continue
      fi


      F0=`grep "Intensive Hamiltonian" $dir3/*out | tail -n 1`
      if [ $? -eq "0" ]; then
        echo -n "$phase " >> ./${filename}
        echo $F0 | awk '{printf(" %s ",$4)}' >> ./${filename}
        echo "" >> ./${filename}
      else
        echo "Warning: no Intensive Hamiltonian in $dir1/$dir2/$dir3/*out "
      fi
    done


    ## Also pull out molecular architecture parameters from the DIS phase
    #if [ ! -e DISPhase/DIS.out ]; then
    #   echo "$dir1/$dir2/DISPhase/DIS.out does not exist!"
    #else
    #    L1=`grep "Arm length" DISPhase/DIS.out | head -n 1 | awk '{print $6}'`
    #    L2=`grep "Arm length" DISPhase/DIS.out | tail -n 1 | awk '{print $6}'`
    #    fA2=`grep "Block fractions" DISPhase/DIS.out | tail -n 1 | awk '{print $6}'`
    #    fB=`grep "Block fractions" DISPhase/DIS.out | tail -n 1 | awk '{print $8}'`
    #    echo "$phiA $L1 $L2 $fA2 $fB" >> ../MiktoParams_vs_phiA.dat
    #fi

    #cd ../
    cd $IDIR/$dir1
  done
  #cd ..
  cd $IDIR
done
