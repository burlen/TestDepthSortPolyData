#!/bin/bash

echo "profiling sort by first point"
valgrind --tool=callgrind --collect-atstart=no  --toggle-collect=vtkDepthSortPolyData2::RequestData'*' --callgrind-out-file=new_by_point.cg ./bin/ds 21 ./data/iso.vtk  ""
echo "profiling sort by bounds centroid"
valgrind --tool=callgrind --collect-atstart=no  --toggle-collect=vtkDepthSortPolyData2::RequestData'*' --callgrind-out-file=new_by_bounds.cg ./bin/ds 22 ./data/iso.vtk  ""
echo "profiling sort by parametric center"
valgrind --tool=callgrind --collect-atstart=no  --toggle-collect=vtkDepthSortPolyData2::RequestData'*' --callgrind-out-file=new_by_pcenter.cg ./bin/ds 23 ./data/iso.vtk  ""

echo "looking for leaks in sort by first point"
valgrind --leak-check=full --track-origins=yes --show-reachable=yes ./bin/ds 21 ./data/iso.vtk  ""
echo "looking for leaks in sort by bounds centroid"
valgrind --leak-check=full --track-origins=yes --show-reachable=yes ./bin/ds 22 ./data/iso.vtk  ""
echo "looking for leaks in sort by parametric center"
valgrind --leak-check=full --track-origins=yes --show-reachable=yes ./bin/ds 23 ./data/iso.vtk  ""
