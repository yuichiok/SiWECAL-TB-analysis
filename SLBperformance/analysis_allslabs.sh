#!/bin/bash

run="Run_ILC_06072020_weekend_3_Ascii"

data_folder="/mnt/win2/Run_Data/"${run}"/"
output="../converter_SLB/convertedfiles/"${run}"/"
mkdir $output

cd ../converter_SLB
root -l -q ConvertDirectorySL.cc\(\"${data_folder}\",false,\"${output}\"\)
hadd ${output}/${run}.root ${output}/*.root
cd -

root -l -q Proto.cc\(\"${output}/${run}\",\"${run}\"\) 


