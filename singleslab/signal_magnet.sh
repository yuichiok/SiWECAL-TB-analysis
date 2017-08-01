#!/bin/bash

dif=dif_1_1_1
path=/home/irles/cernbox/TB2017/TBdata/Magnet/

####################################################################
### REF RUN
#file=${path}0T_ref_3GeV/20170621_143842/run_0_50min_${dif}.raw.root
#output_file=_0T_ref
#root -l SigAnalysisMagnet.cc\(\"$file\",\"$dif\",\"$output_file\"\) 


####################################################################
##1 Tesla Runs
#for i in {0..5}
#do
#    file=${path}1T_3GeV/20170621_164722/run_${i}_60min_${dif}.raw.root
#    output_file=_1T_run$i
#    root -l SigAnalysisMagnet.cc\(\"$file\",\"$dif\",\"$output_file\"\)
#done

j=5
for i in {0..6}
do
    file=${path}1T_3GeV/20170622_005947/run_${i}_60min_${dif}.raw.root
    j=$((j+1))
    output_file=_1T_run${j}
    echo $j
    root -l SigAnalysisMagnet.cc\(\"$file\",\"$dif\",\"$output_file\"\)
done

#######################################################################
##0.5 T runs
for i in {0..3}
do
    file=${path}0.5T_3GeV/20170622_085745/run_${i}_60min_${dif}.raw.root
    output_file=_0.5T_run$i
    root -l SigAnalysisMagnet.cc\(\"$file\",\"$dif\",\"$output_file\"\)
done

######################################################################
#0T post ref
file=${path}0T_post_ref_3GeV/20170622_120758/run_0_60min_${dif}.raw.root
output_file=_0T_post_ref
root -l SigAnalysisMagnet.cc\(\"$file\",\"$dif\",\"$output_file\"\)

