#!/bin/bash

if [ $# -ne 3 ]; then
    echo "Usage: $0 <tau_> <chiN> <threshold to use>"
    exit 1
fi
tau=$1
chiN=$2
threshold=$3

idir=`pwd`

tmpfile="tmp12345"
if [ -e $tmpfile ]; then
    rm $tmpfile
fi


#for i in tau0.9_chiN40.000/phiA0.*/*Phase/bridging/d.bridging/$threshold/; do
for i in *${tau}_chiN${chiN}/phiA0.*/*Phase/bridging/d.bridging/$threshold/; do
    phiA=`echo $i | cut -f 2 -d '/' | sed 's/[a-z,A-Z]//g'`
    phase=`echo $i | cut -f 3 -d '/' `
    
    cd $i
    if [ -e "Pm.log" ]; then
        oneMinusPM=`grep 'Bridging.*1-PM' Pm.log | cut -f 2 -d ':'`
        sumPtilda=`grep 'Bridging.*Ptilda' Pm.log | cut -f 2 -d ':'`
        sumPhat=`grep 'Bridging.*Phat' Pm.log | cut -f 2 -d ':'`
        echo "$phiA $oneMinusPM $sumPtilda $sumPhat $phase" >> $idir/$tmpfile
    #else 
        #echo $i/Pm.log not found

    fi
    cd $idir
done
if [ -e $tmpfile ]; then
    sort -n $tmpfile
    rm $tmpfile
else
    echo "No directories processed..."
fi


