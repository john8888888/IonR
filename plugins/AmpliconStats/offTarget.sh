#!/bin/bash


BAMFILE=$1;
BEDFILE=$2;
OUTPUT_DIR=$3;
BASE_DIR=$4;
REFERENCE=$5;
ANALYSISNAME=$6;
AMPSUM=$7;

OFFTARGETNAME="${OUTPUT_DIR}/${ANALYSISNAME}.offtarget";
OFFTARGETBAM="${OFFTARGETNAME}.bam";
OFFTARGETBED="${OFFTARGETNAME}.bed";
OFFTARGETFA="${OFFTARGETNAME}.fa";
OFFTARGETBEDNAMED="${OFFTARGETNAME}.named.bed";
PROXIMITYBED="${OFFTARGETNAME}.proximity.bed";
OFFTARGETAMPREPORT="${OFFTARGETNAME}.ampSummary.csv"
OFFTARGETSNAP="${OFFTARGETNAME}.snap.tab"
GAPALLOWED=10;
PRIMERLEN=25;
PRIMERINT=10;
TARGETNAME="${OUTPUT_DIR}/${ANALYSISNAME}.target";
TARGETFLANKBED="${TARGETNAME}.flanks.bed";
OFFTARGETFLANKBED="${OFFTARGETNAME}.flanks.bed";
TARGETFLANKFA="${TARGETNAME}.flanks.fa";
OFFTARGETFLANKFA="${OFFTARGETNAME}.flanks.fa";
OFFTARGETALIGN="${OFFTARGETNAME}.flanks.align.sam";
OFFTARGETFLANKBEDALIGNED="${OFFTARGETNAME}.flanks.align.bed";
OFFTARGETFIRST5BASEDAT="${ANALYSISNAME}.offtarget.first.5.prime.base.tab";
ONTARGETFIRST5BASEDAT="${ANALYSISNAME}.ontarget.first.5.prime.base.tab";

echo -n '' > $TARGETFLANKBED;
echo -n '' > $OFFTARGETFLANKBED;
echo -n '' > $OFFTARGETFLANKBEDALIGNED;
echo -n '' > $OFFTARGETBEDNAMED;
echo -n '' > $OFFTARGETSNAP;

echo "[offtargetanalysis] : creating off target bam file";
echo `date`;
echo 'call ion_genOffTargetBam.py';
python $BASE_DIR/ion_genOffTargetBam.py $BAMFILE $BEDFILE | samtools view -Sb - > ${OFFTARGETBAM};
$BASE_DIR/samtools index ${OFFTARGETBAM};
echo "[offtargetanalysis] : creating off target bed file";
echo `date`;
echo 'call ion_buildOffTargetBed.py'; #not long time to run
python $BASE_DIR/ion_buildOffTargetBed.py $OFFTARGETBAM | $BASE_DIR/mergeBed -i - -d $GAPALLOWED > ${OFFTARGETBED};
echo "[offtargetanalysis] : finished creating off target bed file";
echo `date`;

OFFTARGETBEDSIZE=$(stat -c%s "$OFFTARGETBED");
if [ "$OFFTARGETBEDSIZE" != "0" ]; then
    echo "[offtargetanalysis] : creating proximity bed file";
    echo `date`;
    $BASE_DIR/closestBed -d -a ${OFFTARGETBED} -b $BEDFILE > ${PROXIMITYBED}

    #it took an hour to read in bed file of 290k lines, with previous .sh script
    perl $BASE_DIR/ion_offTargetCreateFlankBed.pl $1 $2 $3 $4 $5 $6 $7

    echo `date`;
    echo 'call ion_ampliconReporterOffTarget.py';
    $BASE_DIR/fastaFromBed -name -fi $REFERENCE -bed $OFFTARGETBEDNAMED -fo $OFFTARGETFA;
    OFF='off'
    python $BASE_DIR/ion_ampliconReporterOffTarget.py $OFFTARGETBAM $OFFTARGETBEDNAMED $OFF "no_exist_file" $OFFTARGETAMPREPORT $OFFTARGETFA;

    echo `date`;
    echo 'call ion_generateOffTargetSnapTab.pl';
    perl $BASE_DIR/ion_generateOffTargetSnapTab.pl $OFFTARGETAMPREPORT $OFFTARGETBEDNAMED $OFFTARGETSNAP

fi


#calculate first 5' base
echo -e -n '' > $ONTARGETFIRST5BASEDAT;
echo -e -n '' > $OFFTARGETFIRST5BASEDAT;
$BASE_DIR/count5PrimeBase.sh $BASE_DIR $BAMFILE $REFERENCE >> ${OUTPUT_DIR}/$ONTARGETFIRST5BASEDAT;
$BASE_DIR/count5PrimeBase.sh $BASE_DIR $OFFTARGETBAM $REFERENCE >> ${OUTPUT_DIR}/$OFFTARGETFIRST5BASEDAT;
echo `date`;
echo 'call count5PrimeBase.pl for ontarget';
perl $BASE_DIR/count5PrimeBase.pl $BASE_DIR $BAMFILE $REFERENCE ${OUTPUT_DIR}/$ONTARGETFIRST5BASEDAT;
echo `date`;
echo 'call count5PrimeBase.pl for offtarget';
perl $BASE_DIR/count5PrimeBase.pl $BASE_DIR $OFFTARGETBAM $REFERENCE ${OUTPUT_DIR}/$OFFTARGETFIRST5BASEDAT;

echo `date`;
echo 'call ion_plotOffTargetGC.py';
AMPSUM=$OUTPUT_DIR/$ANALYSISNAME"_ampliconSummary.csv";
python $BASE_DIR/ion_plotOffTargetGC.py $AMPSUM $OFFTARGETAMPREPORT $OUTPUT_DIR