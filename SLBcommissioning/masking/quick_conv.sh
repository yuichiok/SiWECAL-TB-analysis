#!/bin/bash
date=$1
round=$2
it=$3

data_path="/mnt/win2/"
softw_path=$PWD"/../../"

run="Run_ILC_"${date}"_"${round}"_aq"${it}
mkdir ${softw_path}/converter_SLB/convertedfiles/${run}
output=${softw_path}"/converter_SLB/convertedfiles/"${run}"/"

cd ../../converter_SLB
data_folder=${data_path}"/Run_Data/"${run}"/"
root -l -q ConvertDirectorySL_Raw.cc\(\"${data_folder}\",false,\"${run}\",\"${output}\",0\) 
hadd ${output}/../${run}.root ${output}/*.root 
cd -