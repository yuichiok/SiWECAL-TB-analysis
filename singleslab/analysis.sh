#!/bin/bash



for angle in {10,20,30,40}
do
    root -l PedAnalysis.cc\(\"/home/irles/WorkAreaECAL/2018/TB201807/data_converted/LongSlab/angle_sunday8/longslab_angle${angle}_ASU8_dif_1_1_1.raw.root\",\"dif_1_1_1\",\"ASU8_angle${angle}\"\)
    root -l SigAnalysis.cc\(\"/home/irles/WorkAreaECAL/2018/TB201807/data_converted/LongSlab/angle_sunday8/longslab_angle${angle}_ASU8_dif_1_1_1.raw.root\",\"dif_1_1_1\",\"ASU8_angle${angle}\"\)
done

