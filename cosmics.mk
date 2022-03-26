# These variables most likely will have to be changed in your run.
RUN=3GeVMIPscan_run_050043
COMMISSIONING_TAG=3GeVMIPscan_highgain_method2
WOLFRAM_CONFIG=0
SLABS=0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
# COBS=1 11 5 8
COBS=
MAX_ENTRIES_BUILD=-1

# Depending on your setup, you might have to change these as well.
# RAW_DATA_DIR=/eos/project/s/siw-ecal/TB2021-11/commissioning/data/run_050XXX/${RUN}
RAW_DATA_DIR=/eos/project/s/siw-ecal/TB2021-11/beamData/rootfiles/3GeVMIPscan/${RUN}
# COMMISSIONING_FOLDER=/eos/project/s/siw-ecal/TB2021-11/commissioning/results_cosmics
COMMISSIONING_FOLDER=/eos/project/s/siw-ecal/TB2021-11/calibration
# COMMISSIONING_FOLDER=${PWD}

# PEDESTALS=${COMMISSIONING_FOLDER}/pedestals/pedestal_${COMMISSIONING_TAG}.txt
# MIP_CALIB=${COMMISSIONING_FOLDER}/mip_calib/MIP_${COMMISSIONING_TAG}.txt
# MASKED=${COMMISSIONING_FOLDER}/masked/masked_${COMMISSIONING_TAG}.txt

PEDESTALS=${COMMISSIONING_FOLDER}/pedestals/Pedestal_method2_3GeVMIPscan_highgain.txt
MIP_CALIB=${COMMISSIONING_FOLDER}/mip_calib/MIP_3GeVMIPscan_highgain_method2.txt
MASKED=${COMMISSIONING_FOLDER}/masked/masked_channels_3GeVMIPscan_run_050043.txt

# There should be no need to change any of the code below.
CONVERTED_DIR=${PWD}/converter_SLB/convertedfiles/${RUN}
DATA_CONVERTED=${CONVERTED_DIR}/${RUN}_converted.root
TMP=${PWD}/build

build : ${CONVERTED_DIR}/${RUN}_build.root
converted : ${DATA_CONVERTED}
pedestals : ${PEDESTALS}
mip_calib : ${MIP_CALIB}
masked : ${MASKED}

clean:
	rm -rf ${TMP}

# -----------------------------------------------------------------------------
#
# Conversion raw .dat -> .root
#
DAT_PARTS:=$(sort $(patsubst ${RAW_DATA_DIR}/${RUN}.%.tar.gz,%,$(wildcard ${RAW_DATA_DIR}/${RUN}*.tar.gz)))
TMP_CONVERTED_PARTS=${TMP}/converted_parts
TMP_RAW=${TMP}/raw
CONVERTED_PARTS=$(addprefix ${TMP_CONVERTED_PARTS}/converted.,$(addsuffix .root,${DAT_PARTS}))
${TMP_CONVERTED_PARTS}/converted.%.root : ${TMP_RAW}/${RUN}.%
	@mkdir -p ${TMP_CONVERTED_PARTS}
	cd converter_SLB; root -l -q ConvertDataSL.cc\(\"$(word 1,$^)\",false,\"$@\"\)

${TMP_RAW}/${RUN}.% : ${RAW_DATA_DIR}/${RUN}.%.tar.gz
	@mkdir -p ${TMP_RAW}	
	tar xf $(word 1,$^) -C ${TMP_RAW}

${DATA_CONVERTED} : ${CONVERTED_PARTS}
	@echo "[Make] Combine the $(words $^) parts that were converted to .root"
	@mkdir -p $(dir ${DATA_CONVERTED})
	@if [ -f $@ ]; then rm $@; fi
	hadd $@ $(sort $^)

# -----------------------------------------------------------------------------
#
# Conversion _converted.root -> _build.root
#
${CONVERTED_DIR}/${RUN}_build.root : ${DATA_CONVERTED} \
        $(PEDESTALS) $(MIP_CALIB) $(MASKED)
	cd eventbuilding; ./build_events.py $(word 1,$^)\
		--w_config ${WOLFRAM_CONFIG} --cob_positions_string "${COBS}"\
		--out_file_name $@\
		--pedestals_file $(word 2,$^)\
		--mip_calibration_file $(word 3,$^)\
		--masked_file $(word 4,$^)  --max_entries ${MAX_ENTRIES_BUILD}

# Replicate the building procedure, but on each partical converted file -> batch/multiprocessing.
TMP_BUILD_PARTS=${TMP}/build_parts
BUILD_PARTS=$(addprefix ${TMP_BUILD_PARTS}/build.,$(addsuffix .root,${DAT_PARTS}))
${TMP_BUILD_PARTS}/build.%.root : ${TMP_CONVERTED_PARTS}/converted.%.root \
        $(PEDESTALS) $(MIP_CALIB) $(MASKED)
	@mkdir -p ${TMP_BUILD_PARTS}
	cd eventbuilding; ./build_events.py $(word 1,$^)\
		--w_config ${WOLFRAM_CONFIG} --cob_positions_string "${COBS}"\
		--out_file_name $@\
		--pedestals_file $(word 2,$^)\
		--mip_calibration_file $(word 3,$^)\
		--masked_file $(word 4,$^)  --max_entries ${MAX_ENTRIES_BUILD}
		
build_from_parts : ${BUILD_PARTS}
	@if [ ! -f ${CONVERTED_DIR}/${RUN}_build.root ]; then hadd ${CONVERTED_DIR}/${RUN}_build.root $(sort $^); fi

# -----------------------------------------------------------------------------
#
# Read mask from settings
#
${MASKED} : ${TMP_RAW}/Run_Settings.txt
	cd SLBcommissioning; root -l -q -b test_read_masked_channels_summary.C\(\"$(basename $^)\"\)
	mv $(basename $^)_masked.txt $@

${TMP_RAW}/Run_Settings.txt : ${RAW_DATA_DIR}/Run_Settings.txt.tar.gz
	@mkdir -p $(dir $@)
	@tar xf $(word 1,$^) -C $(dir $@)

# -----------------------------------------------------------------------------
#
# Pedestals
#
PEDESTAL_PREFIX=Pedestal_15_layers_
PEDESTAL_ROOT=${PWD}/SLBperformance/results_pedestal/${PEDESTAL_PREFIX}${COMMISSIONING_TAG}.root
# The following lines show the way it is handled in the old shell scripts:
# Build a histogram per partial data file, and then combine them with hadd.
# This hadd step was ridiculously slow, much slower than just rebuilding the hists
# from the full data file.
# This time-tradeoff might change if the data to process becomes much larger.
# Thus, we keep the deprecated code path in the file.
#
# TMP_PEDESTALS=${TMP}/pedestals
# PEDESTAL_PARTS=$(addprefix ${TMP_PEDESTALS}/${PEDESTAL_PREFIX}${COMMISSIONING_TAG}_,$(addsuffix .root,${DAT_PARTS}))
# ${TMP_PEDESTALS}/${PEDESTAL_PREFIX}${COMMISSIONING_TAG}_%.root : ${TMP_CONVERTED_PARTS}/converted.%.root
# 	@echo "[Make] Pedestal analysis preparation with part of the data."
# 	@mkdir -p ${TMP_PEDESTALS}
# 	cd SLBperformance; root -l -q Proto.cc\(\"$(basename $(word 1,$^))\",\"${COMMISSIONING_TAG}_$*\",\"pedestal\"\)
# 	mv SLBperformance/results_proto/${PEDESTAL_PREFIX}${COMMISSIONING_TAG}_$*.root $@
#
# ${PEDESTAL_ROOT} : ${PEDESTAL_PARTS}
# 	hadd $@ $(sort $^)

${PEDESTAL_ROOT} : ${DATA_CONVERTED}
	@echo "[Make] Create pedestal histograms"
	cd SLBperformance; root -l -q Proto.cc\(\"$(basename $(word 1,$^))\",\"${COMMISSIONING_TAG}\",\"pedestal\"\)
	mv SLBperformance/results_proto/${PEDESTAL_PREFIX}${COMMISSIONING_TAG}.root $@

${PEDESTALS} : ${PEDESTAL_ROOT}
	@echo "[Make] Fit pedestals from the histograms"
	cd SLBperformance/results_pedestal; root -l -q PedestalFile.C\(\"$(notdir $(basename $(word 1, $^)))\"\)
	cp SLBperformance/results_pedestal/$(notdir $(basename $(word 1, $^))).txt $@
	@echo "[Make] Done with pedestals"


# -----------------------------------------------------------------------------
#
# MIP calibration
#
MIP_PREFIX=MIPs_15_layers_
MIP_ROOT=${PWD}/SLBperformance/results_mipcalibration/${MIP_PREFIX}${COMMISSIONING_TAG}.root
MIP_ROOT_IN_PROTO_FOLDER=${PWD}/SLBperformance/results_proto/${MIP_PREFIX}${COMMISSIONING_TAG}.root
${MIP_ROOT} : ${DATA_CONVERTED} ${PEDESTALS}
	@echo "[Make] Calibrate MIPs"
	cd SLBperformance; root -l -q Proto.cc\(\"$(basename $(word 1,$^))\",\"${COMMISSIONING_TAG}\",\"mip\",\"$(word 2,$^)\"\)
	cp SLBperformance/results_proto/${MIP_PREFIX}${COMMISSIONING_TAG}.root $@

${MIP_CALIB} : $(addprefix SLBperformance/results_proto/MIPs_layer_,$(addsuffix _${COMMISSIONING_TAG},${SLABS}))
	cp $(word 1,$^) $@
	@# n+3 to skip the header in the subsequent files.
	for layer_calib in $(filter-out $(word 1,$^),$^);do\
		tail -n+3 $${layer_calib} >> $@;\
	done

${MIP_ROOT_IN_PROTO_FOLDER} : ${MIP_ROOT}
	cp ${MIP_ROOT} SLBperformance/results_proto/${MIP_PREFIX}${COMMISSIONING_TAG}.root
	
SLBperformance/results_proto/MIPs_layer_%_${COMMISSIONING_TAG} : ${MIP_ROOT_IN_PROTO_FOLDER}
	@echo "[Make] MIP for slab #$*"
	cd SLBperformance/results_proto; root -l -q -b analysis.cc\(\"${COMMISSIONING_TAG}\",$*\)

# -----------------------------------------------------------------------------
#
# Appendix: Some more scripts as Makefile commands.
#
retriggers : ${DATA_CONVERTED}
	cd SLBperformance; root -l -q Proto.cc\(\"$(basename $(word 1,$^))\",\"${COMMISSIONING_TAG}\",\"retriggers\"\)

monitoring : monitoring_7 monitoring_10

monitoring_% : ${DATA_CONVERTED}
	cd SLBperformance;\
	root -l -q DummyDisplay.cc\(\"$(basename $(word 1,$^))\",\"${COMMISSIONING_TAG}\",$*\)
