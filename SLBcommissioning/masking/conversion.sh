#!/bin/bash
#Expect files in this format:

run="Run_ILC_10152020_3_cosmic_it19_Ascii"

data_folder="/mnt/win2/Run_Data/"${run}"/"

output=${PWD}"/../../converter_SLB/convertedfiles/"${run}"/"
mkdir $output

cd ../../converter_SLB
root -l -q ConvertDirectorySL.cc\(\"${data_folder}\",false,\"${output}\"\)
hadd ${output}/${run}.root ${output}/*.root
cd -

root -l -q TestCheckNoisy.cc
