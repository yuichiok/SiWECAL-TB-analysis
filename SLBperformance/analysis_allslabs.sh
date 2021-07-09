#!/bin/bash

run="run_050001"

data_folder="/eos/project/s/siw-ecal/TB2021-11/commissioning/data/run_050XXX/run_050001_25062021_19h25min_Ascii/"
output="../converter_SLB/convertedfiles/"${run}"/"
mkdir $output

cd ../converter_SLB
root -l -q ConvertDirectorySL.cc\(\"${data_folder}\",false,\"${output}\"\)
hadd ${output}/${run}.root ${output}/*.root
cd -

root -l -q Proto.cc\(\"${output}/${run}\",\"${run}\",\"retriggers\"\) 


