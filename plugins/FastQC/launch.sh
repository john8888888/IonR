#!/bin/bash
# Copyright (C) 2012 Ion Torrent Systems, Inc. All Rights Reserved

VERSION="3.4.1.1" # major.minor.bug

# ===================================================
# Plugin functions
# ===================================================

#*! @function
#  @param  $*  the command to be executed
run ()
{
    echo "running: $*"
    eval "$*";
    EXIT_CODE="$?";
    return ${EXIT_CODE}
}

#*! @function
set_output_paths ()
{
    PLUGIN_OUT_UNMAPPED_BAM_NAME=${TSP_FILEPATH_UNMAPPED_BAM##*/};
}

# ===================================================
# Plugin initialization
# ===================================================
# Test for the existence of the relevant input files.
test_for_file "${TSP_FILEPATH_UNMAPPED_BAM}";

# Set defaults
set_output_paths;

# Test for the existence of the relevant executables.
FASTQC="${DIRNAME}/FastQC/fastqc"
test_for_executable ${FASTQC};

FASTQC_VERSION=`${FASTQC} --version`

# Make sure plugin output is cleared of any previous results
if [ -f ${TSP_FILEPATH_PLUGIN_DIR} ]; then
    rm -rf "${TSP_FILEPATH_PLUGIN_DIR}/*_fastqc"
    rm -rf "${TSP_FILEPATH_PLUGIN_DIR}/*_block.html"
    rm -rf "${TSP_FILEPATH_PLUGIN_DIR}/FastQC_reports.html"
else
    mkdir -p "${TSP_FILEPATH_PLUGIN_DIR}"
fi

# ===================================================
# Run FastQC Plugin
# ===================================================

## Always run on combined file:
FASTQC_OUT_DIR=${PLUGIN_OUT_UNMAPPED_BAM_NAME/%.bam/_fastqc}
run "${FASTQC} ${TSP_FILEPATH_UNMAPPED_BAM} --outdir=${RESULTS_DIR}/"
RET=$?

# Initialize Report
REPORT_HTML=${RESULTS_DIR}/${PLUGINNAME}_block.html
(
    printf "<html>\n<body>\n"
    #Show all data per_base_quality graph on block preview window
    [ ${RET} -eq 0 ] && printf "<img width=\"600\" height=\"450\" src=\"${FASTQC_OUT_DIR}/Images/per_base_quality.png\"/><p>Per Base Sequece Quality - All Reads</p>\n"
    printf "<ul>\n"
    # Add top level report, if successful
    [ ${RET} -eq 0 ] && printf "<li><a target=\"_blank\" href=\"${FASTQC_OUT_DIR}/fastqc_report.html\">FastQC Report - All Reads</a></li>\n"
) > ${REPORT_HTML}

# Barcoded run - append per-run results
if [ -f "${TSP_FILEPATH_BARCODE_TXT}" ]; then
    echo "Found barcode data at: ${TSP_FILEPATH_BARCODE_TXT}"
    NUM_BARCODES=`grep -c barcode ${TSP_FILEPATH_BARCODE_TXT}`;
    echo "${NUM_BARCODES} barcodes found. Report links will appear for barcodes with at least one read in the FASTQ file as analyses complete. You will need to refresh to update.<p/>" >> ${REPORT_HTML};

    for BARCODE_LINE in $(grep "^barcode" ${TSP_FILEPATH_BARCODE_TXT}) "0,nomatch,Non-barcoded" ; do
        BARCODE_ID=$(echo "${BARCODE_LINE}" | cut -d',' -f 2)
        BARCODE_UNMAPPED_BAM="${BARCODE_ID}_${PLUGIN_OUT_UNMAPPED_BAM_NAME}"

        # If that fastq file exists then continue
        if [ -e "${BASECALLER_DIR}/${BARCODE_UNMAPPED_BAM}" ]; then
            run ${FASTQC} "${BASECALLER_DIR}/${BARCODE_UNMAPPED_BAM}" --outdir="${TSP_FILEPATH_PLUGIN_DIR}/";
            RET=$?
            if [ ${RET} -eq 0 ]; then
                # s/.fastq/_fastqc/
                FASTQC_OUT_DIR_BARCODE="${BARCODE_UNMAPPED_BAM%.bam}_fastqc"
                printf "<li><a target=\"_blank\" href=\"${FASTQC_OUT_DIR_BARCODE}/fastqc_report.html\">FastQC Report - ${BARCODE_ID}</a></li>\n" >> ${REPORT_HTML}
            fi
        fi
    done
fi

## Close up Report output
(
    printf "</ul>\n"
    printf "<pre>${FASTQC_VERSION}</pre>\n"
    printf "</body>\n</html>\n"
) >> ${REPORT_HTML}
