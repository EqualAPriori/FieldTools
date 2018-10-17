#!/bin/bash




if [ $# -eq 0 ]; then
    echo "Usage: <list of chiN directories>"
    exit
fi


outfile='bridging_frac.dat'
#for i in `find . -name 'bridging.out'`; do 
#    echo $i; 
#    phiA=`echo $i | cut -f 2 -d '/' | sed 's/[a-z,A-Z]//g'`
#    bridging_frac=`grep 'Bridging fraction:' $i | tail -n 1 | cut -f 2 -d ':'`
#    if [ ! -z $bridging_frac ]; then
#        echo $phiA $bridging_frac >> $outfile
#    fi
#done
idir=`pwd`
#for i in `find . -name 'BridgingOperator.dat'`; do 

for chiNdir in $@; do
    first=1
    cd $chiNdir
    for i in phiA*/*Phase/bridging/d.bridging/BridgingOperator.dat; do 
        echo $chiNdir/$i; 
        mydir=`dirname $i`
        #chiNdir=`echo $i | cut -f 1 -d '/'`
        phi=`echo $i | cut -f 1 -d '/' | sed 's/[a-z,A-Z]//g' | xargs printf %.3f`

        if [ $first -eq 1 ]; then
            rm $idir/$chiNdir/$outfile
            first=0
        fi 

        cd $mydir
        /home/lequieu/miktoarm/SCFT_Mikto/bridging/bridging-scripts/calc_Pm.py > Pm.log
        Ptilda=`grep sumPtilda Pm.log | cut -f 2 -d ' '`
        Phat=`grep sumPhat Pm.log | cut -f 2 -d ' '`
        cd $idir/$chiNdir
        
        if [ -n $Ptilda ] && [ -n $Phat ]; then
            echo $phi $Ptilda $Phat >> $idir/$chiNdir/$outfile
        else
            echo "Error in calc_Pm.py for $mydir, skipping..."
        fi

    done
    cd $idir
done

echo "sorting..."
for i in `find . -name "$outfile"`; do
    echo $i
    sort -n $i > tmp
    mv tmp $i
done
