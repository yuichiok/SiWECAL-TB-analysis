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
 echo "You forgot the 2nd Argument: Round (first1, first2, first3, second, scurves, cosmic)"
 (exit 33)
fi

if [[ $3 == "" ]]
 then
 echo "You forgot the 3rd Argument: iteration"
 (exit 33)
fi


softw_path=$PWD"/../../"
if [ ! -d "${softw_path}/SLBcommissioning/${date}" ]; then
    mkdir ${soft_path}/SLBcommissioning/${date}  
fi
data_path="/mnt/win2/"


if [[ $round == "first1" || $round == "first2" || $round == "first3" ]]
then
    for i in 0 1
    do
	run="Run_ILC_"${date}"_"${round}"_"${i}"_Ascii"
	mkdir ${softw_path}/converter_SLB/convertedfiles/${run}
	output=${softw_path}"/converter_SLB/convertedfiles/"${run}"/"
	
	cd ../../converter_SLB
	data_folder=${data_path}"/Run_Data/"${run}"/"
	root -l -q ConvertDirectorySL.cc\(\"${data_folder}\",false,\"${output}\"\) 
	hadd ${output}/../Run_ILC_${date}_${round}_${i}.root ${output}/*.root 
	cd -
    done

    cp ${data_path}/Run_Data/Run_ILC_${date}_${round}_0_Ascii/Run_Settings.txt  ${softw_path}/SLBcommissioning/${date}/Run_Settings_it${it}.txt
    root -l -q SimpleNoiseChecks.cc\(\"${date}\",\"${round}\",${it}\)
    it2=$((it+1))
    cp ${softw_path}/SLBcommissioning/${date}/Run_Settings_comm_it${it2}.txt ${data_path}/Setup/.
fi


if [[ $round == "second" ]]
then

    run="Run_ILC_"${date}"_"${round}"_"${it}"_Ascii"
    mkdir ${softw_path}"/converter_SLB/convertedfiles/"${run}
    output=${softw_path}"/converter_SLB/convertedfiles/"${run}"/"
    
    cd ../../converter_SLB
    data_folder=${data_path}"/Run_Data/"${run}"/"
    root -l -q ConvertDirectorySL.cc\(\"${data_folder}\",false,\"${output}\"\)
    hadd ${output}/../Run_ILC_${date}_${round}_${it}.root ${output}/*.root
    cd -
    
    cp ${data_path}/Run_Data/Run_ILC_${date}_${round}_${it}_Ascii/Run_Settings.txt  ${softw_path}/SLBcommissioning/${date}/Run_Settings_it${it}.txt
    root -l -q SimpleNoiseChecks.cc\(\"${date}\",\"${round}\",${it}\)
    it2=$((it+1))
    cp ${softw_path}/SLBcommissioning/${date}/Run_Settings_comm_it${it2}.txt ${data_path}/Setup/.
fi

if [[ $round == "scurves" ]]
then
    echo "DID YOU SET UP CORRECTLY THE SLBOARD-ADDress mapping in scurves.C ???"
    cp ${data_path}/Histos/RateVsThresholdScan_${date}_SLBoard.txt ../${date}/.
    cp ${data_path}/Run_Data/Run_ILC_${date}_second_${it}_Ascii/Run_Settings.txt ${softw_path}/SLBcommissioning/${date}/Run_Settings_it${it}.txt
    root -l -q scurves.C\(\"${date}\",${it}\)
    it2=$((it+1))
    cp ${softw_path}/SLBcommissioning/${date}/Run_Settings_comm_it${it2}.txt ${data_path}/Setup/.
    cp histos/scurves_${date}.root ../${date}/.
fi




if [[ $round == "cosmic" ]]
then

    run="Run_ILC_"${date}"_"${round}"_"${it}"_Ascii"
    mkdir ${softw_path}"/converter_SLB/convertedfiles/"${run}
    output=${softw_path}"/converter_SLB/convertedfiles/"${run}"/"

    cd ../../converter_SLB
    data_folder=${data_path}"/Run_Data/"${run}"/"
    root -l -q ConvertDirectorySL.cc\(\"${data_folder}\",false,\"${output}\"\)
    hadd ${output}/../Run_ILC_${date}_${round}_${it}.root ${output}/*.root
    cd -

    cp ${data_path}/Run_Data/Run_ILC_${date}_${round}_${it}_Ascii/Run_Settings.txt  ${softw_path}/SLBcommissioning/${date}/Run_Settings_it${it}.txt
    root -l -q SimpleNoiseChecks.cc\(\"${date}\",\"${round}\",${it}\)
    it2=$((it+1))
    cp ${softw_path}/SLBcommissioning/${date}/Run_Settings_comm_it${it2}.txt ${data_path}/Setup/.
fi
