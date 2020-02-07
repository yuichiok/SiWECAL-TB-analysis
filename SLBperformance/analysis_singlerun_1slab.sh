#!/bin/bash
#Expect files in this format:
# /path/run_XXXXX_SLB_2.root (or SLB_3)
# run example:
#> source analysis.sh XXXXX convert
# where:
#          XXXXX is the run number
#          convert=0 if no conversiom =1 if yes

#source analysis.sh 21010 1

run="Run_ILC_cosmic_test_11212019_11h08min_Ascii"
#Run_ILC_slab14_HV100_12192019_14h58min_Ascii"
data_folder="/home/calice/TB2020/commissioning/data_slboard/"${run}"/"
runname=$run

#cd ../converter_SLB
#root -l -q ConvertDirectorySL.cc\(\"${data_folder}\"\) 
#hadd ${data_folder}/../${runname}.root ${data_folder}/${run}*.root 

#cd -

root -l Monitoring.cc\(\"${data_folder}/../${runname}\",\"${runname}\"\) 

#root -l SingleSlabAnalysis.cc\(\"${data_folder}/../${runname}\",\"${runname}\",0\) 
#root -l SingleSlabAnalysis.cc\(\"${data_folder}/../${runname}\",\"${runname}\",1\) &
#root -l SingleSlabAnalysis.cc\(\"${data_folder}/../${runname}\",\"${runname}\",2\) 

#sleep 10

#cd performance
#source analysis.sh $runname 0
#source analysis.sh $runname 1
#source analysis.sh $runname 2

