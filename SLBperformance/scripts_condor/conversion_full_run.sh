#!/bin/bash

condor=1
datapath="/lustre/ific.uv.es/prj/ific/flc/SiWECAL/TB202206/commissioning/"
runname="cosmic_long_06042022_run_000015"
output=${datapath}"/../converter_SLB/convertedfiles/"${runname}
slbperformancepath=${PWD}"/../"

if [ -d "$output" ]; then                                     
    echo $output 
else                                                          
    mkdir $output                                             
fi 

output2=${PWD}"/../../converter_SLB/convertedfiles/"
if [ -d "$output2" ]; then
    echo $output2
else
    mkdir $output2
fi


for irun in 0 1 2 3 4 5
do
    cat > logs/send_convert_${runname}_${irun}.sh <<EOF
source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-01/init_ilcsoft.sh

data_folder=${datapath}${runname}
echo ${data_folder}
cd ${data_folder}

FILE_start=${runname}"_raw.bin"
if test -f "${FILE_start}"; then
    FILE_new=${FILE_start}"_0000"
    mv ${FILE_start} ${FILE_new}
fi

cd ${slbperformancepath}

cd ../converter_SLB

root -l -q ConvertDirectorySL_Raw.cc\(\"${datapath}/${runname}/\",false,\"${runname}\",\"${datapath}/../converter_SLB/convertedfiles/${runname}\",$irun\) 

cd -

EOF
	    
    cat > logs/send_convert_${runname}_${irun}.sub <<EOF
# Unix submit description file
# send_convert_${runname}_${irun}.sub --
executable              = send_convert_${runname}_${irun}.sh
log                     = send_convert_${runname}_${irun}.log 
output                  = outfile_send_convert_${runname}_${irun}.txt
error                   = errors_send_convert_${runname}_${irun}.txt
should_transfer_files   = Yes
when_to_transfer_output = ON_EXIT
queue 1                          
EOF

    cd logs/
    if (( $condor == 1 )) ; then
	condor_submit send_convert_${runname}_${irun}.sub
	cd -
    else
	source send_convert_${runname}_${irun}.sh &
	cd -
    fi
done
