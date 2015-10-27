#!/bin/bash

echo "timing sort by first point"
tn=`./bin/ds 21 ./iso.vtk  "" 2>&1 | grep = | cut -d= -f2 | cut -d' ' -f2`
to=`./bin/ds 11 ./iso.vtk  "" 2>&1 | grep = | cut -d= -f2 | cut -d' ' -f2`
su=`echo "scale=4; $to/$tn" | bc`
echo "speed up is $su"


echo "timing sort by bounds centroid"
tn=`./bin/ds 22 ./iso.vtk  "" 2>&1 | grep = | cut -d= -f2 | cut -d' ' -f2`
to=`./bin/ds 12 ./iso.vtk  "" 2>&1 | grep = | cut -d= -f2 | cut -d' ' -f2`
su=`echo "scale=4; $to/$tn" | bc`
echo "speed up is $su"


echo "timing sort by parametric center"
tn=`./bin/ds 23 ./iso.vtk  "" 2>&1 | grep = | cut -d= -f2 | cut -d' ' -f2`
to=`./bin/ds 13 ./iso.vtk  "" 2>&1 | grep = | cut -d= -f2 | cut -d' ' -f2`
su=`echo "scale=4; $to/$tn" | bc`
echo "speed up is $su"
