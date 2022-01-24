#!/bin/bash
sleep 1h
for sca in 15 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
do
    for gain in 1 0
    do
        for run in {50043..50124}
	do
	    cat > send_noiseanalysis_${sca}_${gain}_${run}.sh <<EOF
source /lhome/ific/a/airqui/SiWECAL/SiWECAL-TB-analysis-Analysis2021/SLBperformance/send_noise/analysis_run_050xxx_noise.sh 3GeVMIPscan_run_0${run} $gain $sca
EOF
	    
cat > send_noiseanalysis_${sca}_${gain}_${run}.sub <<EOF
# Unix submit description file
# send_noiseanalysis_${sca}_${gain}_${run}.sub --
executable              = send_noiseanalysis_${sca}_${gain}_${run}.sh
log                     = send_noiseanalysis_${sca}_${gain}_${run}.log 
output                  = outfile_send_noiseanalysis_${sca}_${gain}_${run}.txt
error                   = errors_send_noiseanalysis_${sca}_${gain}_${run}.txt
should_transfer_files   = Yes
when_to_transfer_output = ON_EXIT
queue 1                          
EOF

condor_submit send_noiseanalysis_${sca}_${gain}_${run}.sub
	done
    done
done
