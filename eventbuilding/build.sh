for run in {32001..32004}
do
    ./mergeRootFiles.py $run /home/irles/WorkAreaECAL/2019/TB201906/pyrame/run_32/ /home/irles/WorkAreaECAL/2019/TB201906/SLB_data/ /home/irles/WorkAreaECAL/2019/TB201906/EvBuilt/ TDC
    ./build_events.py /home/irles/WorkAreaECAL/2019/TB201906/EvBuilt/run_${run}_merge.root -1 0
done

for run in {32015..32015}
do
    ./mergeRootFiles.py $run /home/irles/WorkAreaECAL/2019/TB201906/pyrame/run_32/ /home/irles/WorkAreaECAL/2019/TB201906/SLB_data/ /home/irles/WorkAreaECAL/2019/TB201906/EvBuilt/ TDC
    ./build_events.py /home/irles/WorkAreaECAL/2019/TB201906/EvBuilt/run_${run}_merge.root -1 0
done


