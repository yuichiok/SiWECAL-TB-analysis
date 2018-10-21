#!/bin/bash

dif="1_1_1"
root -l SigAnalysis.cc\(\"/home/irles/WorkAreaECAL/2018/TB201809/muon_10slab/dif_${dif}.raw.root\",\"dif_$dif\",\"\"\) &
dif="1_1_2"
root -l SigAnalysis.cc\(\"/home/irles/WorkAreaECAL/2018/TB201809/muon_10slab/dif_${dif}.raw.root\",\"dif_$dif\",\"\"\) &
dif="1_1_3"
root -l SigAnalysis.cc\(\"/home/irles/WorkAreaECAL/2018/TB201809/muon_10slab/dif_${dif}.raw.root\",\"dif_$dif\",\"\"\) 
dif="1_1_4"
root -l SigAnalysis.cc\(\"/home/irles/WorkAreaECAL/2018/TB201809/muon_10slab/dif_${dif}.raw.root\",\"dif_$dif\",\"\"\) &
dif="1_1_5"
root -l SigAnalysis.cc\(\"/home/irles/WorkAreaECAL/2018/TB201809/muon_10slab/dif_${dif}.raw.root\",\"dif_$dif\",\"\"\) &

dif="1_2_1"
root -l SigAnalysis.cc\(\"/home/irles/WorkAreaECAL/2018/TB201809/muon_10slab/dif_${dif}.raw.root\",\"dif_$dif\",\"\"\) 
dif="1_2_2"
root -l SigAnalysis.cc\(\"/home/irles/WorkAreaECAL/2018/TB201809/muon_10slab/dif_${dif}.raw.root\",\"dif_$dif\",\"\"\) &
dif="1_2_3"
root -l SigAnalysis.cc\(\"/home/irles/WorkAreaECAL/2018/TB201809/muon_10slab/dif_${dif}.raw.root\",\"dif_$dif\",\"\"\) &
dif="1_2_4"
root -l SigAnalysis.cc\(\"/home/irles/WorkAreaECAL/2018/TB201809/muon_10slab/dif_${dif}.raw.root\",\"dif_$dif\",\"\"\) &
dif="1_2_5"
root -l SigAnalysis.cc\(\"/home/irles/WorkAreaECAL/2018/TB201809/muon_10slab/dif_${dif}.raw.root\",\"dif_$dif\",\"\"\)



