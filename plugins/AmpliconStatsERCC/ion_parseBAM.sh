#!/bin/bash
# Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

if (test $# != 15)
then
    echo "$0 <bam.file> <bed.file> <library> <analysis.name> <amplicon.fasta> <plugindir> <plotSubDir> <hg19ontargetCoverage> <amplicon.summary> <readinfo> <midF0x4SamFile> <biased.base.summary> <merged.bed> <base.coverage.by.chr.tab> <offTargetOn>";
    exit;
fi

BAMFILE=$1;
BEDFILE=$2;
LIBRARY=$3;
ANALYSISNAME=$4;
AMPFASTA=$5;
PLUGINDIR=$6;
PLOTSUBDIR=$7;
#HG19ONTBAM=$8;
#HG19ONTSAM=$8;
#BCOV=$10;
#FCOV=$11;
#RCOV=$12;
DIRECTIONALCOV=$8;
FORBAM=${DIRECTIONALCOV}".for";
REVBAM=${DIRECTIONALCOV}".rev";
AMPSUM=$9;
READINFO=${10};
MIDF0x4Sam=${11};
BIASEDBASE=${12};
MERGEDBED=${13};
CHRSTATS=${14};
OFFTARGETON=${15};
echo "$READINFO"
echo "$MIDF0x4Sam"

#AMPSUM=$PLOTSUBDIR/$ANALYSISNAME"_ampliconSummary.csv";
TARGETSTATS=$PLOTSUBDIR/$ANALYSISNAME".targetStats.txt";
#copy bam file to current location
#cp $BAMFILE copy.bam
echo "current path now is ";
pwd;

echo "generating amplicon summary......"
echo `date`
#generate amplicon summary
ON='on'
#python $PLUGINDIR/ion_ampliconReporter.py $BAMFILE $BEDFILE $ON $MIDF0x4Sam $AMPSUM $AMPFASTA;
BEDSTRING='splitBeds';
if [ -d $BEDSTRING ]; then
    rm -f $BEDSTRING;
else
    mkdir $BEDSTRING;
fi
BEDSTRING=$PLOTSUBDIR$BEDSTRING;
perl $PLUGINDIR/split.bedfile.pl $BEDFILE 10000 50 1000000 $BEDSTRING $BAMFILE
#python $PLUGINDIR/ion_ampliconReporter.py $BAMFILE $BEDFILE $ON $AMPSUM $AMPFASTA;
python $PLUGINDIR/ion_ampliconReporter.py $BAMFILE $BEDFILE $ON $AMPSUM $BEDSTRING $AMPFASTA;
mv split*bed splitBeds;


echo "generating representation images......"
echo `date`
#generate the representation images
#cp $PLUGINDIR/gnuplot.representation.gnuplot $PLOTSUBDIR/;
bedLineNumber=`wc -l $BEDFILE | awk '{print $1}'`;

#PQ='quick';
#if [ $bedLineNumber -gt 1 ]; then
#    echo 'produce quick representation plot'
#    python $PLUGINDIR/ion_plotRepresentation.py $AMPSUM $PLOTSUBDIR $PQ $PLUGINDIR;
#else
#    PQ='high';
#    echo 'produce highqual representation plot'
#    python $PLUGINDIR/ion_plotRepresentation.py $AMPSUM $PLOTSUBDIR $PQ $PLUGINDIR;
#fi
python $PLUGINDIR/ion_plotRepresentation.py $AMPSUM $PLOTSUBDIR $PLUGINDIR;
echo "generating amplicon statistics images......"
echo `date`
#generate the amplicon statistics images
python $PLUGINDIR/ion_plotAmpStats.py $AMPSUM $PLOTSUBDIR;


echo "generating bam stats......"
echo `date`
#get the total number of reads from the run
#if [ -a $BAMFILE".bai" ]; then
#    echo "using index to calculate total reads"
#    $PLUGINDIR/samtools idxstats $BAMFILE | awk 'BEGIN {total=0} {total+=$3; total+=$4 } END {print total}' > $TARGETSTATS;
#else
#    $PLUGINDIR/samtools view $BAMFILE | wc -l > $TARGETSTATS;
#fi
#get the total number of mapped reads
#$PLUGINDIR/samtools view -F 0x4 $BAMFILE | wc -l >> $TARGETSTATS;
#above two stats aleady in readInfo.txt file
head -1 $READINFO > $TARGETSTATS;
#20130415 new requested statistics
#tail -1 $READINFO >> $TARGETSTATS;
#additional line at the end
tail -2 $READINFO | head -1 >> $TARGETSTATS;

#get the total number of ontarget mapped reads
#$PLUGINDIR/samtools view -F 0x4 -L $BEDFILE $BAMFILE -o $PLOTSUBDIR/$HG19ONTSAM;
#wc -l $PLOTSUBDIR/$HG19ONTSAM | cut -f1 -d' ' >> $TARGETSTATS;
#$PLUGINDIR/samtools view -F 0x4 -L $BEDFILE $BAMFILE | wc -l >> $TARGETSTATS; #5m
samtools view -F 0x4 -L $BEDFILE $BAMFILE | wc -l >> $TARGETSTATS; #5m

echo "generating each direction coverage stats......"
echo `date`
#$PLUGINDIR/samtools view -F 16 $BAMFILE -b -o $FORBAM;
#$PLUGINDIR/samtools view -f 16 $BAMFILE -b -o $REVBAM;
#$PLUGINDIR/samtools view -L $BEDFILE -F 16 $BAMFILE -b -o $FORBAM;
#$PLUGINDIR/samtools view -L $BEDFILE -F 20 $BAMFILE -b -o $FORBAM;
#$PLUGINDIR/samtools view -L $BEDFILE -f 16 $BAMFILE -b -o $REVBAM;
#$PLUGINDIR/samtools depth -b $BEDFILE $FORBAM > $FORBAM".depth";
#$PLUGINDIR/samtools depth -b $BEDFILE $REVBAM > $REVBAM".depth";
#version 0.1.19 samtools change -f/-F to -g/-G
#samtools depth -b $BEDFILE $BAMFILE -F 20 > $FORBAM".depth";
#samtools depth -b $BEDFILE $BAMFILE -f 16 > $REVBAM".depth";
samtools depth -b $BEDFILE $BAMFILE -G 20 > $FORBAM".depth";
samtools depth -b $BEDFILE $BAMFILE -g 16 > $REVBAM".depth";
#12m
#awk '{print $1"_"$2"\t"$3}' $FORBAM".depth" | sort > $FORBAM".depth.joined";
#awk '{print $1"_"$2"\t"$3}' $REVBAM".depth" | sort > $REVBAM".depth.joined";
#join -a 1 -a 2 -e '0' -o '0,1.2,2.2' $FORBAM".depth.joined" $REVBAM".depth.joined" | sed -e "s/_/\t/g" | awk '{print $1"\t"$2"\t"$3+$4"\t"$3"\t"$4}' > $DIRECTIONALCOV;
echo "generating both direction stats......"
echo `date`
#$PLUGINDIR/ion_join.coverage.depth.pl $FORBAM".depth" $REVBAM".depth" $DIRECTIONALCOV $BIASEDBASE; #5m
#$PLUGINDIR/ion_join.coverage.depth.pl $FORBAM".depth" $REVBAM".depth" $DIRECTIONALCOV $BIASEDBASE $MERGEDBED $CHRSTATS; #5m
#$PLUGINDIR/ion_join.coverage.depth.pl $FORBAM".depth" $REVBAM".depth" $DIRECTIONALCOV $BIASEDBASE $MERGEDBED $OFFTARGETON $CHRSTATS; #5m
$PLUGINDIR/ion_join.coverage.depth.pl $FORBAM".depth" $REVBAM".depth" $DIRECTIONALCOV $BIASEDBASE $MERGEDBED $OFFTARGETON $CHRSTATS $BEDFILE; #5m
echo `date`


#generate the bam file with ONLY on target reads, for bi-directional coverage analysis
#$PLUGINDIR/samtools view -b -F 0x4 -L $BEDFILE $BAMFILE -o $PLOTSUBDIR/$HG19ONTBAM;
#$PLUGINDIR/samtools view -h $PLOTSUBDIR/$HG19ONTBAM > $PLOTSUBDIR/$HG19ONTSAM;

#generate bi-directional coverage files
#ontarget bam
#$PLUGINDIR/samtools view -F 0x4 -L $BEDFILE $BAMFILE -h -b -o $HG19ONTBAM
#ontarget sam
#$PLUGINDIR/samtools view $HG19ONTBAM -o $HG19ONTSAM

#coverage: both strand, forward strand, reverse strand
#genomeCoverageBed writes temp files to ../../ directory, which is not permitted
#BCOV="both.coverage";
#FCOV="for.coverage";
#RCOV="rev.coverage";
#$PLUGINDIR/genomeCoverageBed -dz -ibam $HG19ONTBAM > $PLOTSUBDIR/$BCOV
#$PLUGINDIR/genomeCoverageBed -dz -strand + -ibam $HG19ONTBAM > $PLOTSUBDIR/$FCOV
#$PLUGINDIR/genomeCoverageBed -dz -strand - -ibam $HG19ONTBAM > $PLOTSUBDIR/$RCOV

