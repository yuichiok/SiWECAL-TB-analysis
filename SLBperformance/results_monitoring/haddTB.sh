#!/bin/bash

name=$1
run=$2

hadd ${name}_${run}.root ${name}_${run}_*.root 

rm ${name}_${run}_*.root 
