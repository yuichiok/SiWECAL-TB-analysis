#!/bin/bash

run=$1
data_folder="/mnt/HardDrive/beamData/ascii/"${run}"/"
output="../converter_SLB/convertedfiles/"${run}"/"

cd $data_folder
for file in *tar.gz
do
    tar xzvf ${file}
done
cd -

initial=${PWD}

cd ${data_folder}
FILE_start=${run}".dat"
if test -f "${FILE_start}"; then
    FILE_new=${FILE_start}"_0000"
    mv ${FILE_start} ${FILE_new}
fi
cd -

if [ -d "$output" ]; then
    cd $output
    ls -1tr | tail -n -1 | xargs -d '\n' rm -f --
    cd -
else
    mkdir $output
fi

cd ../converter_SLB

root -l -q ConvertDirectorySL_TB.cc\(\"${data_folder}\",false,\"${run}\",\"../converter_SLB/convertedfiles/${run}\",2\) &
root -l -q ConvertDirectorySL_TB.cc\(\"${data_folder}\",false,\"${run}\",\"../converter_SLB/convertedfiles/${run}\",1\) &
root -l -q ConvertDirectorySL_TB.cc\(\"${data_folder}\",false,\"${run}\",\"../converter_SLB/convertedfiles/${run}\",0\)

sleep 30s
# source compress_TB2021_endedrun.sh $run &
cd ${initial}

cd ../SLBcommissioning
root -l test_read_masked_channels_TB.C\(\"${data_folder}/Run_Settings.txt\",\"${initial}/../masked/masked_channels_${run}\"\)

rm -rf $data_folder
