# These variables most likely will have to be changed in your run.
RUN=run_050014_07272021_16h35min_Ascii
COMMISSIONING_TAG=_PROTO14_run_050010
WOLFRAM_CONFIG=0

# Depending on your setup, you might have to change these as well.
RAW_DATA_DIR=/eos/project/s/siw-ecal/TB2021-11/commissioning/data/run_050XXX/${RUN}
COMMISSIONING_FOLDER=/eos/project/s/siw-ecal/TB2021-11/commissioning/results_cosmics

# There should be no need to change any of the code below.
CONVERTED_DIR=${PWD}/converter_SLB/convertedfiles/${RUN}
TMP=${PWD}/build
DATA_ZIP_FILES=$(wildcard ${RAW_DATA_DIR}/${RUN}*.tar.gz)

build: ${CONVERTED_DIR}/${RUN}_build.root
converted : ${CONVERTED_DIR}/${RUN}_converted.root

clean:
	rm -rf ${TMP}

TMP_CONVERTED_PARTS=${TMP}/converted_parts
${CONVERTED_DIR}/${RUN}_converted.root : ${DATA_ZIP_FILES}
	@echo "Convert the raw .dat files into a .root tree."
	@echo "The echos from unzipping the $(words $^) .dat files are suppressed."
	@cd ${RAW_DATA_DIR}; for f in $^; do tar xf $$f -C ${RAW_DATA_DIR}; done
	mkdir -p ${TMP}
	mkdir ${TMP_CONVERTED_PARTS}
	cd converter_SLB; root -l -q ConvertDirectorySL.cc\(\"${RAW_DATA_DIR}/\",false,\"${TMP_CONVERTED_PARTS}\"\)
	hadd $@ ${TMP_CONVERTED_PARTS}/* $@

${CONVERTED_DIR}/${RUN}_build.root : ${CONVERTED_DIR}/${RUN}_converted.root
	cd eventbuilding; ./build_events.py $(word 1,$^) ${COMMISSIONING_TAG}\
		-w ${WOLFRAM_CONFIG}\
		-o $@\
		-c ${COMMISSIONING_FOLDER}