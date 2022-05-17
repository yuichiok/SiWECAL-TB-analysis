#!/bin/bash

id=$1

if [[ $1 == "" ]]
then
	echo "You forgot the 1st Argument: run ID --> 050200 - ..."
	(exit 33)
fi

data_path="/eos/project-s/siw-ecal/TB2022-03/beamData/ascii/"

run_name=`find ${data_path} -type d -name "*run_"${id} -printf '%P\n'`

echo "Run name  :" ${run_name}
echo "Full path :" ${data_path}${run_name}

output="./convertedfiles/"${run_name}"/"

if [ ! -d ${output} ] ; then
    mkdir ${output}
fi

data_folder=${data_path}"/"${run_name}"/"
root -l -q ConvertDirectorySL.cc\(\"${data_folder}\",false,\"${output}\"\)
hadd -f ${output}/../${run}.root ${output}/*.root
