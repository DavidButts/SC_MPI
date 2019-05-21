#!/bin/bash

#### Simply a script to load common module for test runs ##########

echo "Purging existing module  and Loading the following Modules:"

purge="module purge"
GCC="module load GCC/7.3.0-2.30"
OPENMPI="module load OpenMPI/3.1.1"
PYTHON="module load Python/3.7.0"

echo $purge
$purge
echo $GCC
$GCC
echo $OPENMPI
$OPENMPI
echo $PYTHON
$PYTHON

module list

