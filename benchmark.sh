#!/bin/bash

echo "+--------------------------------------------+"
echo "| case        | old     | new     | speed up |"
echo "+--------------------------------------------+"

tn=`./bin/ds 21 ./data/iso.vtk  "" 2>&1 | grep = | cut -d= -f2 | cut -d' ' -f2`
to=`./bin/ds 11 ./data/iso.vtk  "" 2>&1 | grep = | cut -d= -f2 | cut -d' ' -f2`
su=`echo "scale=4; $to/$tn" | bc`
echo "| first point | $to | $tn | $su   |"
echo "+--------------------------------------------+"

tn=`./bin/ds 22 ./data/iso.vtk  "" 2>&1 | grep = | cut -d= -f2 | cut -d' ' -f2`
to=`./bin/ds 12 ./data/iso.vtk  "" 2>&1 | grep = | cut -d= -f2 | cut -d' ' -f2`
su=`echo "scale=4; $to/$tn" | bc`
echo "| axes bounds | $to | $tn | $su   |"
echo "+--------------------------------------------+"

tn=`./bin/ds 23 ./data/iso.vtk  "" 2>&1 | grep = | cut -d= -f2 | cut -d' ' -f2`
to=`./bin/ds 13 ./data/iso.vtk  "" 2>&1 | grep = | cut -d= -f2 | cut -d' ' -f2`
su=`echo "scale=4; $to/$tn" | bc`
echo "| param. cen. | $to | $tn | $su   |"
echo "+--------------------------------------------+"
