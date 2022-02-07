#!/bin/bash
#Expect files in this format:

# run="Run_ILC_10152020_3_cosmic_it19_Ascii"
# run="Run_ILC_16062021_cosmiclong_it17_Ascii"
# run="run_050001_25062021_19h25min_Ascii"
run="Run_ILC_01232022_masking_it0_0_Ascii"

data_folder="/mnt/win2/Run_Data/"${run}"/"

output=${PWD}"/../../converter_SLB/convertedfiles/"${run}"/"
mkdir $output

cd ../../converter_SLB
root -l -q ConvertDirectorySL.cc\(\"${data_folder}\",false,\"${output}\"\)
hadd ${output}/${run}.root ${output}/*.root
cd -

# root -l -q TestCheckNoisy.cc
root -l -q TestCheckNoisy.cc\(\"${output}/${run}.root\",\"${data_folder}/Run_Settings.txt\",\"${output}/Run_Settings_comm_it19.txt\"\)
