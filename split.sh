#!/bin/bash
echo -n "Splitting data..."
split -b 50M iso.vtk data/iso.vtk.
echo "done!"
