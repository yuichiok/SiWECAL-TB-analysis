#!/bin/bash
condor=1
datapath="/lustre/ific.uv.es/prj/ific/flc/SiWECAL/TB202206/commissioning/"
slbperfomancepath=${PWD}"/../"
runname="cosmic_long_06042022_run_000015"
scriptspath=${PWD}


for i in 0 1 2 3 4 5 6 7
do

    cd ${scriptspath}/
    cat > logs/send_mippedestal_${runname}_${i}.sh <<EOF
#!/bin/bash 
    
source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-01/init_ilcsoft.sh
cd ${slbperfomancepath}

hadd -f -k ${datapath}/../converter_SLB/convertedfiles/${runname}_merged_${i}.root ${datapath}/../converter_SLB/convertedfiles/${runname}/*${i}?.root 

cd ${slbperfomancepath}
root -l -q Proto.cc\(\"${datapath}/../converter_SLB/convertedfiles/${runname}_merged_${i}\",\"${runname}_merged_${i}\",0\) 

#cd TBchecks
#source analysis.sh ${runname}_merged_${i}


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
