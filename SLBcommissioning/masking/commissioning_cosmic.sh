#!/bin/bash

date=$1
round=$2
it=$3
nslabs=$4

if [[ $1 == "" ]]
 then
 echo "You forgot the 1st Argument: date --> 05182020 for the 5th May 2020"
 (exit 33)
fi

if [[ $2 == "" ]]
 then
 echo "You forgot the 2nd Argument: Round (masking, scurves, cosmic)"
 (exit 33)
fi

if [[ $3 == "" ]]
 then
 echo "You forgot the 3rd Argument: iteration"
 (exit 33)
fi

softw_path=$PWD"/../../"
if [ ! -d "${softw_path}/SLBcommissioning/${date}" ]; then
    mkdir ${softw_path}/SLBcommissioning/${date}  
fi
if [ ! -d "${softw_path}/SLBcommissioning/masking/histos" ]; then
    mkdir ${softw_path}/SLBcommissioning/masking/histos
fi

# data_path="/mnt/win1/"
data_path="/mnt/win2/"


if [[ $round == "cosmic" ]]
then

    run="Run_ILC_"${date}"_"${round}"_it"${it}"_Ascii"
    mkdir ${softw_path}"/converter_SLB/convertedfiles/"${run}
    output=${softw_path}"/converter_SLB/convertedfiles/"${run}"/"

    cd ../../converter_SLB
    data_folder=${data_path}"/Run_Data/"${run}"/"
    root -l -q ConvertDirectorySL.cc\(\"${data_folder}\",false,\"${output}\"\)
    hadd ${output}/../Run_ILC_${date}_${round}_it${it}.root ${output}/*.root
    cd -

    cp ${data_path}/Run_Data/Run_ILC_${date}_${round}_it${it}_Ascii/Run_Settings.txt  ${softw_path}/SLBcommissioning/${date}/Run_Settings_it${it}.txt
    runfile=${softw_path}/SLBcommissioning/${date}/Run_Settings_it${it}.txt
    it2=$((it+1))
    runfile_new=${softw_path}/SLBcommissioning/${date}/Run_Settings_comm_it${it2}.txt
    root -l -q TestCheckNoisy.cc\(\"${output}/../Run_ILC_${date}_${round}_it${it}.root\",\"${runfile}\",\"${runfile_new}\"\)
    cp ${softw_path}/SLBcommissioning/${date}/Run_Settings_comm_it${it2}.txt ${data_path}/Setup/.
fi
