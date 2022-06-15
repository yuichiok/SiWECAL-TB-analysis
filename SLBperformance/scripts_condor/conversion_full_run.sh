#!/bin/bash

runN=$1
condor=0
datapath="/home/calice/TB2022-06/eudaq_data/raw/"
runname="raw_siwecal_"${runN}
output=${PWD}"/../../converter_SLB/convertedfiles/"${runname}
slbperformancepath=${PWD}"/../"
start=$PWD

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


for irun in 0
do
    cat > logs/send_convert_${runname}_${irun}.sh <<EOF

data_folder=${datapath}${runname}

cd ${slbperformancepath}

cd ../converter_SLB
root -l -q ConvertDirectorySL_Raw_TB202206.cc\(\"${datapath}/\",false,\"${runname}\",\"${output2}/${runname}\",$irun\) 

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
	source send_convert_${runname}_${irun}.sh 
	#cd $start
	#source mippedestal_full_run.sh ${runN} 0 &
    fi
done
cd $start
