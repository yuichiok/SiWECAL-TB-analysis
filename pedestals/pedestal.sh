#!/bin/bash
## For all grid points, merge all runs with same grid number and change naming to have two numbers (00-80)

dif=1_1_1
root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/dif_$dif.raw.root\",\"dif_$dif\",\"\"\) &
dif=1_1_2
root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/dif_$dif.raw.root\",\"dif_$dif\",\"\"\) &
dif=1_1_3
root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/dif_$dif.raw.root\",\"dif_$dif\",\"\"\) &
dif=1_1_4
root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/dif_$dif.raw.root\",\"dif_$dif\",\"\"\) 

dif=1_1_5
root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/dif_$dif.raw.root\",\"dif_$dif\",\"\"\) &
dif=1_2_1
root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/dif_$dif.raw.root\",\"dif_$dif\",\"\"\) &
dif=1_2_2
root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/dif_$dif.raw.root\",\"dif_$dif\",\"\"\) 


for grid in {0..80}
do
    if [ $grid -lt 10 ]
    then
	dif=1_1_1
	root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/grid0${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid0${grid}\"\) &
	dif=1_1_2
	root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/grid0${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid0${grid}\"\) &
	dif=1_1_3
	root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/grid0${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid0${grid}\"\) &
	dif=1_1_4
	root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/grid0${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid0${grid}\"\) 
	dif=1_1_5
	root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/grid0${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid0${grid}\"\) &
	dif=1_2_1
	root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/grid0${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid0${grid}\"\) &
	dif=1_2_2
	root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/grid0${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid0${grid}\"\) 
    else
	dif=1_1_1
	root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/grid${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid${grid}\"\) &
	dif=1_1_2
	root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/grid${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid${grid}\"\) &
	dif=1_1_3
	root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/grid${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid${grid}\"\) &
	dif=1_1_4
	root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/grid${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid${grid}\"\) 
	dif=1_1_5
	root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/grid${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid${grid}\"\) &
	dif=1_2_1
	root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/grid${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid${grid}\"\) &
	dif=1_2_2
	root -l PedAnalysis.cc\(\"/home/irles/WorkArea/TB2017/TBdata/MIPscan/rootfiles/grid${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid${grid}\"\) 
    fi
done
