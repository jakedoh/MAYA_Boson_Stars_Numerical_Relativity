#!/bin/bash

XORIGIN=-50
YORIGIN=-50
XMAX=50
YMAX=50
Z=0
NX=101
NY=101
NMAX=10
LMAX=4
MMAX=4

HDF5FILE="/tmp/psi4_Decomp.h5"
name=$(basename ${HDF5FILE})
name=${name%%.h5}

ITERATION_START=0
ITERATION_END=10000
DELTA_ITERATION=1000

for (( iteration=${ITERATION_START};\
       iteration<= ${ITERATION_END};
       iteration+=${DELTA_ITERATION}))
do
  ./h5read ${HDF5FILE} ${iteration} \
         ${XORIGIN} ${YORIGIN} ${XMAX} ${YMAX} ${Z} ${NX} ${NY} \
         ${NMAX} ${LMAX} ${MMAX} > ${name}.${iteration}.asc
done

