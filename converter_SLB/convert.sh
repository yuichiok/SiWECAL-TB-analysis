#!/bin/bash

date=$1
round=$2
it=$3

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

# data_path="/mnt/win1/"
data_path="/mnt/win2/"


if [[ $round == "cosmic" ]]
then
	run="Run_ILC_"${date}"_"${round}"_it"${it}"_Ascii"
	mkdir "./convertedfiles/"${run}
	output="./convertedfiles/"${run}"/"

	run="Run_ILC_"${date}"_"${round}"_it"${it}"_Ascii"
	data_folder=${data_path}"/Run_Data/"${run}"/"
	root -l -q ConvertDirectorySL.cc\(\"${data_folder}\",false,\"${output}\"\)
	hadd -f ${output}/../Run_ILC_${date}_${round}_it${it}.root ${output}/*.root
fi

if [[ $round == "pedestal" ]]
then
	run="Run_ILC_"${date}"_"${round}"_ch"${it}"_Ascii"
	mkdir "./convertedfiles/"${run}
	output="./convertedfiles/"${run}"/"

	run="Run_ILC_"${date}"_"${round}"_ch"${it}"_Ascii"
	data_folder=${data_path}"/Run_Data/"${run}"/"
	root -l -q ConvertDirectorySL.cc\(\"${data_folder}\",false,\"${output}\"\)
	hadd -f ${output}/../Run_ILC_${date}_${round}_ch${it}.root ${output}/*.root
fi