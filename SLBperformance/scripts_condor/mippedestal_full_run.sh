#!/bin/bash

runN=$1
condor=0
datapath="/home/calice/TB2022-06/eudaq_data/raw/"
runname="raw_siwecal_"${runN}
slbperfomancepath=${PWD}"/../"
scriptspath=${PWD}


for i in 0
do

    cd ${scriptspath}/
    cat > logs/send_mippedestal_${runname}_${i}.sh <<EOF
#!/bin/bash 
    
cd ${slbperfomancepath}
root -l -q Proto.cc\(\"../converter_SLB/convertedfiles/${runname}/converted_${runname}.raw\",\"${runname}\",0\) 

#cd TBchecks
#source analysis.sh ${runname}


EOF

    cd ${scriptspath}/	    
    cat > logs/send_mippedestal_${runname}_${i}.sub <<EOF
# Unix submit description file
# send_mippedestal_${runname}_${i}.sub --
executable              = send_mippedestal_${runname}_${i}.sh
log                     = send_mippedestal_${runname}_${i}.log 
output                  = outfile_send_mippedestal_${runname}_${i}.txt
error                   = errors_send_mippedestal_${runname}_${i}.txt
should_transfer_files   = Yes
when_to_transfer_output = ON_EXIT
queue 1                          
EOF

    echo "submit"
    cd ${scriptspath}/logs/
    if (( $condor == 1 )) ; then
	condor_submit send_mippedestal_${runname}_${i}.sub
	cd -
    else
	source send_mippedestal_${runname}_${i}.sh &
	cd -
    fi
done
