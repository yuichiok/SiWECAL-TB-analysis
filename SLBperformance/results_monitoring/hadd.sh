#!/bin/bash

name=$1
run=$2

hadd ${name}_${run}_0.root ${name}_${run}_0*.root &
hadd ${name}_${run}_1.root ${name}_${run}_1*.root 
hadd ${name}_${run}_2.root ${name}_${run}_2*.root

hadd ${name}_${run}.root ${name}_${run}_0.root ${name}_${run}_1.root ${name}_${run}_2.root

rm ${name}_${run}_0.root ${name}_${run}_1.root ${name}_${run}_2.root
