#!/bin/bash

run=$1

for layer in {0..14}
do
    hadd Retriggers_Layer_${layer}_${run}_2.root Retriggers_Layer_${layer}_*_${run}_2*.root
    hadd Retriggers_Layer_${layer}_${run}_1.root Retriggers_Layer_${layer}_*_${run}_1*.root
    hadd Retriggers_Layer_${layer}_${run}_0.root Retriggers_Layer_${layer}_*_${run}_0*.root

    hadd Retriggers_Layer_${layer}_${run}.root Retriggers_Layer_${layer}_${run}_0.root Retriggers_Layer_${layer}_${run}_1.root Retriggers_Layer_${layer}_${run}_2.root
    rm Retriggers_Layer_${layer}_${run}_0.root Retriggers_Layer_${layer}_${run}_1.root Retriggers_Layer_${layer}_${run}_2.root
done

