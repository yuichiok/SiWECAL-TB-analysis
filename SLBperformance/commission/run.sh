#!/bin/bash

# run example
# run = "Run_ILC_03012022_cosmic_it17"

run=$1
output=${PWD}"/../../converter_SLB/convertedfiles/"

cd ../
root -l -q Proto.cc\(\"${output}${run}\",\"${run}\"\)
cd -