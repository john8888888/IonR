#!/bin/bash

PLUGINDIR=$1;
BAMFILE=$2;
#REFERENCEFASTA="/results/referenceLibrary/tmap-f3/hg19/hg19.fasta"
REFERENCEFASTA=$3;
FIRSTBASEBEDFILE="first.5.prime.base.bed"
GRANDT=0;
TTOTAL=0;
ATOTAL=0;
CTOTAL=0;
GTOTAL=0;


function showCounts {
    DIRECTION=$1
    $PLUGINDIR/ion_calcFirstBaseBed.py ${BAMFILE} ${DIRECTION} > ${FIRSTBASEBEDFILE}
    
    #samtools view ${BAMFILE} | awk '{if ($2 == 0){print $3"\t"$4-2"\t"$4-1"\t"$1}}' > ${FIRSTBASEBEDFILE}
    TOTALREADS=`$PLUGINDIR/samtools view -${DIRECTION} 0x10 ${BAMFILE} | wc -l`;
    TCOUNTS=`$PLUGINDIR/fastaFromBed -s -name -tab -fi ${REFERENCEFASTA} -bed ${FIRSTBASEBEDFILE} -fo stdout | cut -f2 | grep T | wc -l`
    ACOUNTS=`$PLUGINDIR/fastaFromBed -s -name -tab -fi ${REFERENCEFASTA} -bed ${FIRSTBASEBEDFILE} -fo stdout | cut -f2 | grep A | wc -l`
    CCOUNTS=`$PLUGINDIR/fastaFromBed -s -name -tab -fi ${REFERENCEFASTA} -bed ${FIRSTBASEBEDFILE} -fo stdout | cut -f2 | grep C | wc -l`
    GCOUNTS=`$PLUGINDIR/fastaFromBed -s -name -tab -fi ${REFERENCEFASTA} -bed ${FIRSTBASEBEDFILE} -fo stdout | cut -f2 | grep G | wc -l`
    let TTOTAL=$TTOTAL+$TCOUNTS;
    let ATOTAL=$ATOTAL+$ACOUNTS;
    let CTOTAL=$CTOTAL+$CCOUNTS;
    let GTOTAL=$GTOTAL+$GCOUNTS;
    if [ "$DIRECTION" == "f" ]; then
	echo -e "Reverse\t$TCOUNTS\t$ACOUNTS\t$CCOUNTS\t$GCOUNTS";
    else
	echo -e "Forward\t$TCOUNTS\t$ACOUNTS\t$CCOUNTS\t$GCOUNTS";
    fi
}
echo -e "Direction\tT\tA\tC\tG";
#echo "Forward Read Stats"
showCounts f
#echo "Reverse Read Stats"
showCounts F
echo -e "Total\t$TTOTAL\t$ATOTAL\t$CTOTAL\t$GTOTAL";
let GRANDT=$TTOTAL+$ATOTAL+$CTOTAL+$GTOTAL;
TPERC=`echo -e "$TTOTAL $GRANDT" | awk '{print $1*100/$2}'`;
APERC=`echo -e "$ATOTAL $GRANDT" | awk '{print $1*100/$2}'`;
CPERC=`echo -e "$CTOTAL $GRANDT" | awk '{print $1*100/$2}'`;
GPERC=`echo -e "$GTOTAL $GRANDT" | awk '{print $1*100/$2}'`;

echo -e "Percent\t$TPERC\t$APERC\t$CPERC\t$GPERC";
