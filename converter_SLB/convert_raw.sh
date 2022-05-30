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
data_path="/eos/project/s/siw-ecal/TB2022-06/commissioning/"
run="Run_ILC_"${date}"_"${round}"_it"${it}
data_folder=${data_path}"/"${run}"/"

cd ${data_folder}
FILE_start=${run}"_raw.bin"
if test -f "${FILE_start}"; then
FILE_new=${FILE_start}"_0000"
mv ${FILE_start} ${FILE_new}
fi

cd -

output="./convertedfiles/"${run}
if [ -d "$output" ]; then
cd $output
ls -1tr | tail -n -1 | xargs -d '\n' rm -f --
cd -
else
mkdir $output
fi

echo $data_folder

if [[ $round == "cosmic" ]]
then
	root -l -q ConvertDirectorySL_Raw.cc\(\"${data_folder}\",false,\"${run}\",\"${output}\",0\)
	hadd -f ${output}/../Run_ILC_${date}_${round}_it${it}.root ${output}/*.root
fi

if [[ $round == "pedestal" ]]
then
	root -l -q ConvertDirectorySL_Raw.cc\(\"${data_folder}\",false,\"${run}\",\"${output}\",0\)
	hadd -f ${output}/../Run_ILC_${date}_${round}_ch${it}.root ${output}/*.root
fi

if [[ $round == "run" ]]
then
	root -l -q ConvertDirectorySL_Raw.cc\(\"${data_folder}\",false,\"${run}\",\"${output}\",0\)
	hadd -f ${output}/../${run}.root ${output}/*.root
fi