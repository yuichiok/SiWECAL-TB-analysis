#!/bin/bash
## For all grid points, merge all runs with same grid number and change naming to have two numbers (00-80)

sleep 5h

for bcid in {5,10,15,30,60}
do
    for grid in {0..81}
    do
	if [ $grid -lt 10 ]
	then
	    dif=1_1_1
	    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/grid0${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid0${grid}_bcidTh${bcid}\"\) &
	    dif=1_1_2
	    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/grid0${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid0${grid}_bcidTh${bcid}\"\) &
	    dif=1_1_3
	    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/grid0${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid0${grid}_bcidTh${bcid}\"\) &
	    dif=1_1_4
	    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/grid0${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid0${grid}_bcidTh${bcid}\"\) 
	    dif=1_1_5
	    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/grid0${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid0${grid}_bcidTh${bcid}\"\) &
	    dif=1_2_1
	    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/grid0${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid0${grid}_bcidTh${bcid}\"\) &
	    dif=1_2_2
	    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/grid0${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid0${grid}_bcidTh${bcid}\"\) 
	else
	    dif=1_1_1
	    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/grid${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid${grid}_bcidTh${bcid}\"\) &
	    dif=1_1_2
	    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/grid${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid${grid}_bcidTh${bcid}\"\) &
	    dif=1_1_3
	    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/grid${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid${grid}_bcidTh${bcid}\"\) &
	    dif=1_1_4
	    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/grid${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid${grid}_bcidTh${bcid}\"\) 
	    dif=1_1_5
	    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/grid${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid${grid}_bcidTh${bcid}\"\) &
	    dif=1_2_1
	    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/grid${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid${grid}_bcidTh${bcid}\"\) &
	    dif=1_2_2
	    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/grid${grid}_dif_$dif.raw.root\",\"dif_$dif\",\"grid${grid}_bcidTh${bcid}\"\) 
	fi
    done
    
    dif=1_1_1
    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/dif_$dif.raw.root\",\"dif_$dif\",\"bcidTh${bcid}\"\) &
    dif=1_1_2
    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/dif_$dif.raw.root\",\"dif_$dif\",\"bcidTh${bcid}\"\) &
    dif=1_1_3
    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/dif_$dif.raw.root\",\"dif_$dif\",\"bcidTh${bcid}\"\) &
    dif=1_1_4
    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/dif_$dif.raw.root\",\"dif_$dif\",\"bcidTh${bcid}\"\) 
    
    dif=1_1_5
    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/dif_$dif.raw.root\",\"dif_$dif\",\"bcidTh${bcid}\"\) &
    dif=1_2_1
    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/dif_$dif.raw.root\",\"dif_$dif\",\"bcidTh${bcid}\"\) &
    dif=1_2_2
    root -l PedAnalysis.cc\(\"/home/irles/cernbox/TB2017/TBdata/MIPscan/rootfiles_bcidTh${bcid}/dif_$dif.raw.root\",\"dif_$dif\",\"bcidTh${bcid}\"\) 
done

