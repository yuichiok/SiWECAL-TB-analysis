#!/bin/bash
run=$1
run_file="converted"
#output="/eos/project/s/siw-ecal/TB2021-11/beamData/rootfiles/3GeVMIPscan/"${run}"/"
output=${PWD}"/../converter_SLB/convertedfiles/"${run}"/"


initial_folder=$PWD

cd $initial_folder

hadd -f ${output}${run}_merged.root ${output}converted_${run}.d*.root
root -l -q Proto.cc\(\"${output}${run}_merged\",\"${run}\"\)  

