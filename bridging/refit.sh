#!/bin/bash

for i in tau0.9_chiN40.000/phiA0.*/*Phase/bridging/d.bridging/; do 
    cd $i
    filename=threshold_vs_bfraction.dat 
    if [ ! -e $filename ]; then
        echo $i/$filename does not exist!
        continue
    fi


    awk '{if (NF > 1) print $1}' $filename > good_thresholds.dat
    begin=`head -n1 good_thresholds.dat`
    end=`tail -n1 good_thresholds.dat`

gnuplot << EOF
f(x)=a*x+b
fit [$begin:$end] f(x) './$filename' u 1:2 via a,b
set print "bridging_fraction.dat"
print b
EOF

done
