#!/bin/bash
#sleep 3h
source pedestal.sh
cd macros_pedestal/
root -l -q Analysis_bcid.C\(\"dif_1_1_1\"\)
root -l -q Analysis_bcid.C\(\"dif_1_1_2\"\)
root -l -q Analysis_bcid.C\(\"dif_1_1_3\"\)
root -l -q Analysis_bcid.C\(\"dif_1_1_4\"\)
root -l -q Analysis_bcid.C\(\"dif_1_1_5\"\)
root -l -q Analysis_bcid.C\(\"dif_1_2_1\"\)
root -l -q Analysis_bcid.C\(\"dif_1_2_2\"\)
cd ../
source signal.sh
