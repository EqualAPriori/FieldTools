#!/bin/bash

j=0;
for i in `ls density*bin | sed 's/[a-z,A-Z,.]//g' | sort -n`; do 
    echo $i, $j
    jpad=`printf %08g $j`
    vtkfile=density_${jpad}.vtk
    if [ ! -e $vtkfile ]; then
        PolyFTS_to_VTK.py density${i}.bin; 
        mv density${i}.vtk $vtkfile; 
    else
        echo "$vtkfile already exists. Skipping..."
    fi
    let j+=1; 
done

