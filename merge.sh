#!/bin/bash

echo -n "Unpackaging data..."
cd data
cat iso.vtk.aa  iso.vtk.ab  iso.vtk.ac  iso.vtk.ad  iso.vtk.ae > ../iso.vtk
cd ..
echo "done!"
