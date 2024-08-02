#!/bin/bash
# Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

ulimit -s 8192

VERSION="0.4.4"

# grab PUI parameters
PLUGIN_TARGETS=$PLUGINCONFIG__UNPTARGETS
echo "PLUGIN_TARGETS is $PLUGIN_TARGETS";


PLUGINCONFIG__ANALYSISTYPE='Resume'; #default
PLUGINCONFIG__OFFTARGET='Off'; #default

#if [ -z "$PLUGINCONFIG__UNPTARGETS" ]; then
    OLD_IFS="$IFS"
    IFS=";"
#    PLAN_INFO=(`${DIRNAME}/parse_plan.py ${TSP_FILEPATH_PLUGIN_DIR}/startplugin.json.self`);
    PLAN_INFO=(`${DIRNAME}/parse_plan.py ${TSP_FILEPATH_PLUGIN_DIR}/startplugin.json`);
    PLUGINCONFIG__LIBRARYTYPE=${PLAN_INFO[0]};
    #PLUGINCONFIG__VARIATIONTYPE=${PLAN_INFO[1]};
    PLUGINCONFIG__TARGETREGIONS=${PLAN_INFO[1]}; #bed file
    #PLUGINCONFIG__TARGETLOCI=${PLAN_INFO[3]};
    PLUGINCONFIG__ANALYSISTYPE=${PLAN_INFO[2]}; #analysis type
    PLUGINCONFIG__OFFTARGET=${PLAN_INFO[3]}; #off target analysis
    PLUGINCONFIG__PREVNUMBER=${PLAN_INFO[4]}; #previous same-name plugin output number

    echo "library type is $PLUGINCONFIG__LIBRARYTYPE";
    echo "analysis type is $PLUGINCONFIG__ANALYSISTYPE";
    PLUGINCONFIG__LIBRARYTYPE__LC=`echo $PLUGINCONFIG__LIBRARYTYPE | awk '{print tolower($0)}'`;
    echo "previous same-name plugin output number is $PLUGINCONFIG__PREVNUMBER";

if [ -z "$PLUGINCONFIG__UNPTARGETS" ]; then
    if [ -z "$PLUGINCONFIG__TARGETREGIONS" ]; then
        rm -f "${TSP_FILEPATH_PLUGIN_DIR}/results.json"
        HTML="${TSP_FILEPATH_PLUGIN_DIR}/${PLUGINNAME}.html"
        echo '<html><body>' > "$HTML"
        echo "<h3><center>${PLUGIN_RUN_NAME}</center></h3>" >> "$HTML"
        echo "<br/><h2 style=\"text-align:center;color:red\">*** A bed file is required to launch the plugin automatically. ***</h2>" >> "$HTML"
        echo "<br/><h3 style=\"text-align:center\">(Requires an associated Plan that is not a GENS Runtype.)</h3></br>" >> "$HTML"
        echo '</body></html>' >> "$HTML"
        exit
    else
	#if [[ "$PLUGINCONFIG__LIBRARYTYPE__LC" =~ ampliseq ]]; then
        #PLUGIN_TARGETS=`echo "$PLUGINCONFIG__TARGETREGIONS" | sed -e 's/\/unmerged\/detail//'`;
        PLUGIN_TARGETS=`echo "$PLUGINCONFIG__TARGETREGIONS"`;
	echo "Plugin launched automatically via run plan."
	#fi;
    fi
else
    #get it from PUI
    PLUGIN_TARGETS=$PLUGINCONFIG__UNPTARGETS
    echo "Plugin launched manually thru PUI."
fi

#get previous plugin output number
#echo "api url: ${RUNINFO__API_URL}"
#echo "plugin name: ${RUNINFO__PLUGIN_NAME}"
#`${DIRNAME}/ion_get.previous.plugin.output.number.pl ${RUNINFO__API_URL} ${RUNINFO__PLUGIN_NAME}`
#exit;

# Check for merged BAM file
if [ -n "$PLUGINCONFIG__MERGEDBAM" ]; then
    TSP_FILEPATH_BAM=$PLUGINCONFIG__MERGEDBAM
fi

TIMESTAMP=`date +%s`;
if [ -f 'SummaryTable.txt' ]
then 
    mv SummaryTable.txt SummaryTable.txt.${TIMESTAMP}
fi
if [ -f 'results.json' ]
then 
    mv results.json results.json.${TIMESTAMP}
fi
if [ -f 'AmpliconStats.html' ]
then 
    mv AmpliconStats.html AmpliconStats.html.${TIMESTAMP}
fi


if [ ! -f $TSP_FILEPATH_BARCODE_TXT ]
then
    echo "analysis-dir $ANALYSIS_DIR";
    echo "out-dir $RESULTS_DIR";
    echo "pre-out $PLUGINCONFIG__ANALYSISTYPE";
    echo "off-target $PLUGINCONFIG__OFFTARGET";
    echo "bed-file $PLUGIN_TARGETS";
    echo "bam-file $TSP_FILEPATH_BAM";
    echo "web-bam $TSP_URLPATH_BAM";
    echo "fastq $TSP_FILEPATH_FASTQ";
    echo "analysis-name $TSP_ANALYSIS_NAME";
    echo "library $TSP_LIBRARY";
    echo "run-name $TSP_RUN_NAME";
    echo "version $VERSION";
    echo "[ampgen: processing a non-barcode run]";
    PRE_OUT='y'

                $DIRNAME/AmpliconStats.pl \
                        --analysis-dir $ANALYSIS_DIR \
                        --out-dir $RESULTS_DIR \
                        --pre-out $PLUGINCONFIG__ANALYSISTYPE \
                        --pre-num $PLUGINCONFIG__PREVNUMBER \
                        --off-target $PLUGINCONFIG__OFFTARGET \
		        --bed-file $PLUGIN_TARGETS \
		        --bam-file $TSP_FILEPATH_BAM \
		        --reference $TSP_FILEPATH_GENOME_FASTA \
		        --web-bam $TSP_URLPATH_BAM \
			--fastq $TSP_FILEPATH_FASTQ \
			--analysis-name $TSP_ANALYSIS_NAME \
			--library $TSP_LIBRARY \
			--run-name $TSP_RUN_NAME \
                        --version $VERSION
else
    echo "a barcoded run, calling barcode.sh";
        $DIRNAME/barcode.sh $PLUGIN_TARGETS $VERSION $PLUGINCONFIG__ANALYSISTYPE $PLUGINCONFIG__OFFTARGET $PLUGINCONFIG__PREVNUMBER
fi
