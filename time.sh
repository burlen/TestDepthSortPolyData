#!/bin/bash

function time_case {
code=$1
meth=$2
let code_idx=$1-1
let meth_idx=$2-1
code_names=('old' 'new')
meth_names=('first_point' 'cell_bound_cen' 'param_cen')

t='('
t+=`./bin/ds ${code}${meth} ./data/iso.vtk  "" 2>&1 | grep = | cut -d= -f2 | cut -d' ' -f2`
t+=' + '
t+=`./bin/ds ${code}${meth} ./data/iso.vtk  "" 2>&1 | grep = | cut -d= -f2 | cut -d' ' -f2`
t+=' + '
t+=`./bin/ds ${code}${meth} ./data/iso.vtk  "" 2>&1 | grep = | cut -d= -f2 | cut -d' ' -f2`
t+=')/3.0'

echo ${code_names[$code_idx]}_${meth_names[$meth_idx]}=$t
}

time_case $1 1
time_case $1 2
time_case $1 3
