#!/bin/bash
# Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved
# $Revision: 17949 $


#PLUGIN_TARGETS=$PLUGINCONFIG__UNPTARGETS
#echo "PLUGIN_TARGETS is $PLUGIN_TARGETS";
PLUGIN_TARGETS=$1;
VERSION=$2;
PLUGINCONFIG__ANALYSISTYPE=$3
PLUGINCONFIG__OFFTARGET=$4
PLUGINCONFIG__PREVNUMBER=$5
#VERSION="0.2.8"

echo "[AmpStats: processing a barcode run]";
echo "[AmpStats: running on `hostname`"];

echo -n "[analysis type: ]";
echo $PLUGINCONFIG__ANALYSISTYPE
echo -n "[off-target option: ]";
echo $PLUGINCONFIG__OFFTARGET

BARCODE_ZIP_FILE="$ANALYSIS_DIR/*.barcode.bam.zip"
BARCODE_BAI_ZIP_FILE="$ANALYSIS_DIR/*.barcode.bam.bai.zip"
BARCODE_FASTQ_FILE="$ANALYSIS_DIR/*.barcode.fastq.zip"
BARCODE_DIR="$RESULTS_DIR";
BARCODE_WORK_DIR="$BARCODE_DIR/to-delete";
HTML_NAME="${RESULTS_DIR}/AmpliconStats.html";
JSON_NAME="${RESULTS_DIR}/results.json";
HTML_SUMMARY_NAME="${RESULTS_DIR}/SummaryTable.txt";
FIELD_FILE="${DIRNAME}/summary_field.txt";
REFERENCE=$TSP_FILEPATH_GENOME_FASTA;

EXPORT_CSV=$TSP_RUN_NAME"_AmpliconStats.csv";
EXPORT_TAB="AmpliconStats.tab";
echo -n '' > $EXPORT_TAB;


ln -sf ${DIRNAME}/js ${TSP_FILEPATH_PLUGIN_DIR}/.;
ln -sf ${DIRNAME}/css ${TSP_FILEPATH_PLUGIN_DIR}/.;
ln -sf ${DIRNAME}/export ${TSP_FILEPATH_PLUGIN_DIR}/.;

writeUnsuitableHtml(){
	echo -n '' > ${HTML_NAME};
	writeHtmlHeader >> ${HTML_NAME};
	echo "<h1>This run is not suitable for this plugin.</h1>" >> ${HTML_NAME}
	writeHtmlFooter >> ${HTML_NAME};
}

#add export csv script p1
writeHtmlHeader() {
	echo \
	"<html>
	<head>
		<title>AmpliconStats</title>
                <script type="text/javascript" src="/site_media/jquery/js/jquery-1.7.1.min.js"></script>
                <script type="text/javascript" src="./js/sorttable.js"></script>

                <link rel="stylesheet" type="text/css" href="/site_media/stylesheet.css" />
                <link type="text/css" href="/site_media/jquery/css/aristo/jquery-ui-1.8.7.custom.css" rel="stylesheet" />
                <link href="/site_media/jquery/js/tipTipX/jquery.tipTipX.css" rel="stylesheet" type="text/css" />
                <link rel="stylesheet" type="text/css" href="/site_media/report.css" media="screen" />"

        cat ${DIRNAME}"/export_csv_html_p1.txt"
	echo \
	"       </head>
	<body>
                <style type=\"text/css\">
                        table {!important;border-collapse:collapse;margin:auto;table-layout:fixed}
                        tr.nobarcode td {color:#A0A0A0}
                        th,td {font-family:\"Lucida Sans Unicode\",\"Lucida Grande\",Sans-Serif;font-size:14px;line-height:1.2em;font-weight:normal}
                        th,td {!important;border:1px solid #bbbbbb;padding:5px;text-align:center}
                        td {padding-top:20px;padding-bottom:20px}
                        img.frm {display:block;margin-left:auto;margin-right:auto;margin-top:10px;margin-bottom:10px;width:400;height:400;border-width:0;cursor:help}
                        .thelp {cursor:help}
                </style>

		<h1><center>Amplicon General Analysis with AmpliconStats Plugin Version ${VERSION}<br>${TSP_ANALYSIS_NAME}</center></h1>";
}

#add export csv script p2
writeHtmlFooter() {
    echo -n \
"
<br/>
<form id='exportform' method='POST' action='"
    echo -n $TSP_URLPATH_PLUGIN_DIR
    echo -n \
"/export/export.php'>
<input id='exportdata' name='exportdata' type='hidden'/>
<input id='exportfn' name='exportfn' type='hidden' value='";
    echo -n $EXPORT_CSV
    echo \
"'/>
</form>";
    cat ${DIRNAME}"/export_csv_html_p2.txt"
    echo \
	"	</body>
	 </html>";
}

writeTableHeader() {
    echo "<tr>";
    if [ -f "$FIELD_FILE" ]; then
	echo "<th>Barcode</th>";
	while read line
	do
	    echo "<th>$line</th>";
	done < $FIELD_FILE
    fi
    echo "</tr>"
}

writeColumnHeader() {
    if [ -f "$FIELD_FILE" ]; then
	echo -e -n "Barcode";
	while read line
	do
	    echo -e -n "\t$line";
	done < $FIELD_FILE
    fi
}

writeHTML(){
	echo $1;
	echo -n '' > ${HTML_NAME};
	writeHtmlHeader >> ${HTML_NAME};

	#echo -n "<button type='button' id='export'>Export To CSV</button>" >> ${HTML_NAME};
	echo  "<br><a href=\"$EXPORT_TAB\">Download</a><br>" >> ${HTML_NAME}

	echo " Click column headers to sort <br>"  >> ${HTML_NAME};
	echo "<table  id='metrics_tab' name='metrics_tab' class='sortable' border=1 cellpadding=6 width=100%>" >> ${HTML_NAME};

	writeTableHeader >> ${HTML_NAME};

	cat ${HTML_SUMMARY_NAME} >> ${HTML_NAME};
	echo "</table>">> ${HTML_NAME};
	if [ ! -z "$1" ]; then
	echo \
	"
		<h1><center> <font color=#FF0000> ANALYSIS HAS NOT COMPLETED</font></center></h1>
		<center> <font color=#000000> $1</font></center>
		<center> <font color=#0000FF> (Hit Refresh to Update Results)</font></center>
	" >> ${HTML_NAME};
	fi
	writeHtmlFooter >> ${HTML_NAME};

}

#writeHTMLTop;

echo -n '' > ${HTML_SUMMARY_NAME};
#echo "<tr>" >> ${HTML_SUMMARY_NAME};
#if [ -f "$FIELD_FILE" ]; then
#    echo "<th>Barcode</th>" >> ${HTML_SUMMARY_NAME};
#    while read line
#    do
#	echo "<th>$line</th>" >> ${HTML_SUMMARY_NAME};
#    done < $FIELD_FILE
#fi
#echo "</tr>" >> ${HTML_SUMMARY_NAME}

writeHTML "(Copying files...please wait)";	

if [ ! -d $BARCODE_WORK_DIR ]; then 
	mkdir -p $BARCODE_WORK_DIR; 
fi
# now make sure we have all the .bam and .fastq files
cd $BARCODE_WORK_DIR #to-delete
if [ -f $BARCODE_ZIP_FILE ]; then
	#echo "unzipping bam files";
	unzip -n $BARCODE_ZIP_FILE #REMEMBER to un-comment this line
	if [ -f $BARCODE_BAI_ZIP_FILE ]; then
	    #echo "unzipping bai files";
	    unzip -n $BARCODE_BAI_ZIP_FILE #REMEMBER to un-comment this line
	fi;
else
	#echo "copying across all .bam files to temporary storage";
	#cp ${ANALYSIS_DIR}/*.bam . #REMEMBER to un-comment this line
	#cp ${ANALYSIS_DIR}/*.bam.bai . #REMEMBER to un-comment this line
	for i in ${ANALYSIS_DIR}/*.bam ${ANALYSIS_DIR}/*.bam.bai; do
	    ln -s $i
	done
fi

cd $BARCODE_DIR #results_dir, i.e. plugin_out/AmpliconStats_out
for barcode in `grep ^barcode $ANALYSIS_DIR/barcodeList.txt|cut -f2 -d ','`;do
        echo "working on $barcode directory"
	if [ ! -d "$barcode" ]; then 
		mkdir -p "$barcode";
	else
	        if [ -f $barcode"/just.html.txt" ]; then
		    echo $barcode"/just.html.txt";
		else
		    rm -r $barcode/*AmpStats*
		fi
	fi

	#BARCODE_BAM_FILE=$BARCODE_WORK_DIR/$barcode"_R_*.bam";
	#BARCODE_BAI_FILE=$BARCODE_WORK_DIR/$barcode"_R_*.bam.bai";
	BARCODE_BAM_FILE=$BARCODE_WORK_DIR/$barcode"_*.bam";
	BARCODE_BAI_FILE=$BARCODE_WORK_DIR/$barcode"_*.bam.bai";
	if [ -f $BARCODE_BAM_FILE ]; then
	        cd $barcode
		#mv $BARCODE_BAM_FILE "$barcode/"
		#mv $BARCODE_BAI_FILE "$barcode/"
		ln -s $BARCODE_WORK_DIR/$barcode"_rawlib.bam"
		ln -s $BARCODE_WORK_DIR/$barcode"_rawlib.bam.bai"
		#echo "$BARCODE_BAM_FILE linked to $barcode/"
		cd ..
	fi
done

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


writeColumnHeader >> $EXPORT_TAB;



#echo "checking if special barcode/bedfile specified ...";
BARCODEBEDFILEPAIR=`grep -A10 pluginconfig ${TSP_FILEPATH_PLUGIN_DIR}/startplugin.json | grep barcodetargetregions | cut -f2 -d: | cut -f2 -d'"'`;
FIRSTBEDINPAIR=`echo $BARCODEBEDFILEPAIR | cut -f1 -d';' | cut -f2 -d=`;
if [ "${#FIRSTBEDINPAIR}" -gt "5" ]; then BBSPECIFIED=true; else BBSPECIFIED=false; fi;
echo "barcode bedfile pair: $BARCODEBEDFILEPAIR so BBSPECIFIED=$BBSPECIFIED";

echo -e -n "{\n  \"barcoded\":\"true\",\n  \"bedfile\":\"" > ${JSON_NAME};
echo -e $PLUGIN_TARGETS"\",\n  \"barcodes\":{\n" >> ${JSON_NAME}
for barcode in `grep ^barcode $ANALYSIS_DIR/barcodeList.txt|cut -f2 -d ',' | sort`; do
#	BAMNAME=$BARCODE_DIR/$barcode/$barcode"_R_*.bam";
#	BAMNAME=$BARCODE_DIR/$barcode/$barcode"*.bam";
	BAMNAME=$BARCODE_DIR/$barcode/$barcode"*rawlib.bam";
	OUTDIR="$RESULTS_DIR/$barcode";
	writeHTML "(Analysing barcode: $barcode...please wait)";
	echo "[AmpStats: running analysis on $barcode]";
	echo "<tr>">> ${HTML_SUMMARY_NAME};

#	if [[ "$barcode" == IonXpress_00* || "$barcode" == IonXpress_011 || "$barcode" == IonXpress_012 || "$barcode" == IonXpress_013 ]]; then
#	    continue;
#	fi


	if ($BBSPECIFIED); then
	    CURRENTBEDFILE=`echo $BARCODEBEDFILEPAIR | awk -F $barcode '{print $1"\t"$2}' | cut -f2 | cut -f1 -d';' | cut -f2 -d'='`;
	    #echo "current bed file should be $CURRENTBEDFILE";
	    if [ "${#CURRENTBEDFILE}" -gt "5" ]; then 
		HASTHEBARCODE=true; 
	    else 
		HASTHEBARCODE=false; 
	    fi

	    if ($HASTHEBARCODE); then
		PLUGIN_TARGETS=$CURRENTBEDFILE;
	    else
		continue;
	    fi
	#else
	    #echo "no barcode bedfile specification"
	    #no extra action
	fi

	echo "HASTHEBARCODE is $HASTHEBARCODE, bed file changed to $CURRENTBEDFILE";



	#if [ -e $BAMNAME ]; then #REMEMBER to un-comment this line

            echo $ANALYSIS_DIR;
            echo $RESULTS_DIR;
            echo $PLUGIN_TARGETS;
	    echo "reference is $REFERENCE";
	    echo "library is $TSP_LIBRARY";
	    echo "pre-out $PLUGINCONFIG__ANALYSISTYPE";
	    echo "off-target $PLUGINCONFIG__OFFTARGET";
            echo "$ANALYSIS_DIR/*.bam";
            echo "$ANALYSIS_DIR/*.fastq";
            echo $VERSION
	    PRE_OUT='y'


            $DIRNAME/AmpliconStats.pl \
                --analysis-dir $ANALYSIS_DIR \
                --out-dir $OUTDIR \
                --pre-out $PLUGINCONFIG__ANALYSISTYPE \
                --pre-num $PLUGINCONFIG__PREVNUMBER \
                --off-target $PLUGINCONFIG__OFFTARGET \
                --bed-file $PLUGIN_TARGETS \
                --bam-file $BAMNAME \
		--reference $REFERENCE \
		--library $TSP_LIBRARY \
		--web-bam $BAMNAME \
		--analysis-name $barcode \
                --version $VERSION

	    RETVAL=1
	    if [ -f $BARCODE_DIR/$barcode/return.value.txt ]; then
		RETVAL=`cat $BARCODE_DIR/$barcode/return.value.txt`
	    fi
            if (test $RETVAL = 0); then
                    echo "<td><a href="${barcode}/AmpliconStats.html">$barcode</a></td>">> ${HTML_SUMMARY_NAME};
		    echo -e -n "\n$barcode" >> $EXPORT_TAB;
                    if [ -f "${barcode}/summary_table.txt" ]; then
                        summary=($(cat ${barcode}/summary_table.txt | cut -f2))
                        for i in "${summary[@]}"
                        do
                            echo "<td>$i</td>" >> ${HTML_SUMMARY_NAME}
			    echo -e -n "\t$i" >> $EXPORT_TAB;
                        done
		    fi
		    echo -e -n "    \"$barcode\":" >> ${JSON_NAME};
                    if [ -f "${barcode}/results.json" ]; then
                        jsonLine=($(cat ${barcode}/results.json | cut -f2))
                        for i in "${jsonLine[@]}"
                        do
			    if [ "$i" == "}" ]; then
				echo -n "      $i" >> ${JSON_NAME};
			    else
				echo "      $i" >> ${JSON_NAME}
			    fi
                        done
		    fi
		    echo "," >> ${JSON_NAME};
            fi

        #else
                #echo "<td>$barcode<br>(no bam file)</td>">> ${HTML_SUMMARY_NAME};
	#fi #REMEMBER to un-comment this line
	echo "</tr>">> ${HTML_SUMMARY_NAME};
done
echo -e "  }\n}" >> ${JSON_NAME};

if [ -f ${JSON_NAME} ]; then
    JSON_LINE=`wc -l ${JSON_NAME} | cut -f1 -d' '`;
    let TOP_PART=$JSON_LINE-3;
    head -${TOP_PART} ${JSON_NAME} > ${JSON_NAME}".temp";
    tail -3 ${JSON_NAME} | sed -e 's/\(.*\),/\1/' >> ${JSON_NAME}".temp";
    mv ${JSON_NAME}".temp" ${JSON_NAME};
fi



writeHTML;
#rm -rf $BARCODE_WORK_DIR

#generate statistics of whole run 
perl $DIRNAME/ion_generateWholeRunStatistics.pl $BARCODE_DIR
