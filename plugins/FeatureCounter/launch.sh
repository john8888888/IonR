#!/bin/bash
# Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved

# Change the following line to all CAPS to disable auto-run of this plugin, but do not uncomment
#autorundisable

ulimit -s 8192

export PYTHONPATH=$PYTHONPATH:${DIRNAME}/HTSeq-0.5.3p3/lib/python

VERSION="1.0.6" # major.minor.bug


# ===================================================
# Plugin functions
# ===================================================

#*! @function
#  @param  $*  the command to be executed
run ()
{
        echo "running: $*" >> ${TSP_FILEPATH_PLUGIN_DIR}/${PLUGIN_OUT_DETAILED_LOG};
        eval $*;
        EXIT_CODE="$?";
}

#*! @function
set_output_paths ()
{
        PLUGIN_OUT_BAM_NAME=`echo ${TSP_FILEPATH_BAM} | sed -e 's_.*/__g'`; 
        PLUGIN_OUT_DETAILED_LOG=${PLUGINNAME}_log.txt;
}

# set defaults
set_output_paths;

# Make sure it is empty
if [ -f ${TSP_FILEPATH_PLUGIN_DIR} ]; then
    run "rm -rf ${TSP_FILEPATH_PLUGIN_DIR}";
fi

#BED_FILE=$PLUGINCONFIG__unptargets;
echo "$PLUGINCONFIG__UNPTARGETS";
#MODE=$PLUGINCONFIG__mode;
echo "$PLUGINCONFIG__MODE";

perl $DIRNAME/scripts/bed2gff.pl $PLUGINCONFIG__UNPTARGETS $TSP_FILEPATH_PLUGIN_DIR/flag.txt > $TSP_FILEPATH_PLUGIN_DIR/features.gff
out=`cat $TSP_FILEPATH_PLUGIN_DIR/flag.txt`;
echo "$out";

#barcoded run
if [ -f ${TSP_FILEPATH_BARCODE_TXT} ]; then
    echo "${TSP_FILEPATH_BARCODE_TXT}";

    BARCODE_LINES=`wc -l ${TSP_FILEPATH_BARCODE_TXT} | awk '{print \$1}'`;
    echo ${BARCODE_LINES};
    NUM_BARCODES=`expr ${BARCODE_LINES} - 2`;

    echo "Report links will appear as analyses complete on each barcode." > ${TSP_FILEPATH_PLUGIN_DIR}/samples_block.html;
    
    CTR=0;
    for BARCODE_LINE in `cat ${TSP_FILEPATH_BARCODE_TXT} | grep "^barcode"`
    do
        BARCODE_ID=`echo ${BARCODE_LINE} | awk 'BEGIN{FS=","} {print $2}'`;
        BARCODE_SEQ=`echo ${BARCODE_LINE} | awk 'BEGIN{FS=","} {print $3}'`;
        BARCODE_BAM_NAME="${BARCODE_ID}_${PLUGIN_OUT_BAM_NAME}";
	
	if [ -f ${ANALYSIS_DIR}/${BARCODE_BAM_NAME} ]; then
        	#create sub dir
        	if [ ! -f ${TSP_FILEPATH_PLUGIN_DIR}/${BARCODE_ID}.${BARCODE_SEQ} ]; then
            	run "mkdir -p ${TSP_FILEPATH_PLUGIN_DIR}/${BARCODE_ID}.${BARCODE_SEQ}";
        	fi

        	#create sym link to bam file
        	if [ ! -f ${TSP_FILEPATH_PLUGIN_DIR}/${BARCODE_ID}.${BARCODE_SEQ}/${BARCODE_BAM_NAME} ]; then
            	run "ln -s ${ANALYSIS_DIR}/${BARCODE_BAM_NAME} ${TSP_FILEPATH_PLUGIN_DIR}/${BARCODE_ID}.${BARCODE_SEQ}/${BARCODE_BAM_NAME}";
	    	run "ln -s ${TSP_FILEPATH_PLUGIN_DIR}/features.gff ${TSP_FILEPATH_PLUGIN_DIR}/${BARCODE_ID}.${BARCODE_SEQ}/features.gff";
        	fi


		echo "samtools view -h ${TSP_FILEPATH_PLUGIN_DIR}/${BARCODE_ID}.${BARCODE_SEQ}/${BARCODE_BAM_NAME} | $DIRNAME/HTSeq-0.5.3p3/scripts/htseq-count -q $out -m $PLUGINCONFIG__MODE - $TSP_FILEPATH_PLUGIN_DIR/${BARCODE_ID}.${BARCODE_SEQ}/features.gff > $TSP_FILEPATH_PLUGIN_DIR/${BARCODE_ID}.${BARCODE_SEQ}/stats.txt";
	
		run "samtools view -h ${TSP_FILEPATH_PLUGIN_DIR}/${BARCODE_ID}.${BARCODE_SEQ}/${BARCODE_BAM_NAME} | $DIRNAME/HTSeq-0.5.3p3/scripts/htseq-count -q $out -m $PLUGINCONFIG__MODE - $TSP_FILEPATH_PLUGIN_DIR/${BARCODE_ID}.${BARCODE_SEQ}/features.gff > $TSP_FILEPATH_PLUGIN_DIR/${BARCODE_ID}.${BARCODE_SEQ}/stats.txt";

		echo "$DIRNAME/scripts/grab_stats.pl $TSP_FILEPATH_PLUGIN_DIR/${BARCODE_ID}.${BARCODE_SEQ}/stats.txt $TSP_FILEPATH_PLUGIN_DIR/${BARCODE_ID}.${BARCODE_SEQ}/features.gff $TSP_FILEPATH_PLUGIN_DIR/${BARCODE_ID}.${BARCODE_SEQ}/feature_counts.txt ${TSP_FILEPATH_PLUGIN_DIR}/${BARCODE_ID}.${BARCODE_SEQ}/${BARCODE_BAM_NAME} $TSP_FILEPATH_PLUGIN_DIR/${BARCODE_ID}.${BARCODE_SEQ} > $TSP_FILEPATH_PLUGIN_DIR/${BARCODE_ID}.${BARCODE_SEQ}/summary_stats.txt";

		run "perl $DIRNAME/scripts/grab_stats.pl $TSP_FILEPATH_PLUGIN_DIR/${BARCODE_ID}.${BARCODE_SEQ}/stats.txt $TSP_FILEPATH_PLUGIN_DIR/${BARCODE_ID}.${BARCODE_SEQ}/features.gff $TSP_FILEPATH_PLUGIN_DIR/${BARCODE_ID}.${BARCODE_SEQ}/feature_counts.txt ${TSP_FILEPATH_PLUGIN_DIR}/${BARCODE_ID}.${BARCODE_SEQ}/${BARCODE_BAM_NAME} $TSP_FILEPATH_PLUGIN_DIR/${BARCODE_ID}.${BARCODE_SEQ} > $TSP_FILEPATH_PLUGIN_DIR/${BARCODE_ID}.${BARCODE_SEQ}/summary_stats.txt";
	
		run "perl $DIRNAME/scripts/write_html.pl $TSP_FILEPATH_PLUGIN_DIR ${DIRNAME} ${URL_ROOT} ${BARCODE_ID}.${BARCODE_SEQ}";

        	CTR=`expr ${CTR} + 1`;
     	fi
    done


#non barcoderun

else

if [ ! -f ${TSP_FILEPATH_PLUGIN_DIR}/FeatureCounter ]; then
        run "mkdir -p ${TSP_FILEPATH_PLUGIN_DIR}/FeatureCounter";
fi

if [ ! -f ${TSP_FILEPATH_PLUGIN_DIR}/FeatureCounter/${PLUGIN_OUT_BAM_NAME} ]; then
        run "ln -s ${TSP_FILEPATH_BAM} ${TSP_FILEPATH_PLUGIN_DIR}/FeatureCounter/${PLUGIN_OUT_BAM_NAME}";
	run "ln -s ${TSP_FILEPATH_PLUGIN_DIR}/features.gff ${TSP_FILEPATH_PLUGIN_DIR}/FeatureCounter/features.gff";
fi

echo "samtools view -h $TSP_FILEPATH_BAM | $DIRNAME/HTSeq-0.5.3p3/scripts/htseq-count -q $out -m $PLUGINCONFIG__MODE - $TSP_FILEPATH_PLUGIN_DIR/FeatureCounter/features.gff > $TSP_FILEPATH_PLUGIN_DIR/FeatureCounter/stats.txt";

	run "samtools view -h $TSP_FILEPATH_BAM | $DIRNAME/HTSeq-0.5.3p3/scripts/htseq-count -q $out -m $PLUGINCONFIG__MODE - $TSP_FILEPATH_PLUGIN_DIR/features.gff > $TSP_FILEPATH_PLUGIN_DIR/FeatureCounter/stats.txt";

echo "$DIRNAME/scripts/grab_stats.pl $TSP_FILEPATH_PLUGIN_DIR/FeatureCounter/stats.txt $TSP_FILEPATH_PLUGIN_DIR/features.gff $TSP_FILEPATH_PLUGIN_DIR/FeatureCounter/feature_counts.txt ${TSP_FILEPATH_PLUGIN_DIR}/FeatureCounter/${PLUGIN_OUT_BAM_NAME} $TSP_FILEPATH_PLUGIN_DIR/FeatureCounter > $TSP_FILEPATH_PLUGIN_DIR/FeatureCounter/summary_stats.txt";

	run "perl $DIRNAME/scripts/grab_stats.pl $TSP_FILEPATH_PLUGIN_DIR/FeatureCounter/stats.txt $TSP_FILEPATH_PLUGIN_DIR/FeatureCounter/features.gff $TSP_FILEPATH_PLUGIN_DIR/FeatureCounter/feature_counts.txt ${TSP_FILEPATH_PLUGIN_DIR}/FeatureCounter/${PLUGIN_OUT_BAM_NAME} $TSP_FILEPATH_PLUGIN_DIR/FeatureCounter > $TSP_FILEPATH_PLUGIN_DIR/FeatureCounter/summary_stats.txt";
	
	run "perl $DIRNAME/scripts/write_html.pl $TSP_FILEPATH_PLUGIN_DIR ${DIRNAME} ${URL_ROOT} FeatureCounter";

fi
