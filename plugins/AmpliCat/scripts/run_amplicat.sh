#!/bin/bash
# Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved

#--------- Begin command arg parsing ---------

CMD=`echo $0 | sed -e 's/^.*\///'`
DESCR="Wrapper for ampliCat scripts to also handle read trimming and report workups.
Most of the options passed are not implented but are reserved for future usage.
Individual plots and data files are produced to the output directory ('.' unless specified by -D).
An HTML file for visualizing all results in a browser is also produced, unless the -x option is used."
USAGE="USAGE:
 $CMD [options] <reference.fasta> <BAM file>"
OPTIONS="OPTIONS:
  -h --help Report usage and help
  -a Customize output for Amplicon reads.
  -d Filter to remove Duplicate reads removed.
  -u Filter to Uniquely mapped reads (SAM MAPQ>0).
  -i Apply to Isolated (non-overlapping) targets only.
  -r Customize output for AmpliSeq-RNA reads. (Overides -a.)
  -t Filter BAM file to trimmed reads using TRIMP.
  -p <int>  Padding value for BED file padding. For reporting only. Default: 0.
  -A <file> Annotate coverage for (annotated) targets specified in this BED file
  -B <file> Limit coverage to targets specified in this BED file
  -C <name> Original name for BED targets selected for reporting (pre-padding, etc.)
  -D <dirpath> Path to root Directory where results are written. Default: ./
  -F <file> Create this Filtered BAM file (on-target reads that are not considered as due to false priming. Default: ''.
  -G <file> Genome file. Assumed to be <reference.fasta>.fai if not specified.
  -O <file> Output file name for text data (per analysis). Default: '' => <BAMROOT>.stats.txt.
  -P <file> Padded targets BED file for padded target coverage analysis
  -Q <file> Name for BLOCK HTML results file (in output directory). Default: '' (=> none created)
  -R <file> Name for HTML Results file (in output directory). Default: 'results.html'
  -S <file> SampleID tracking regions file. Default: '' (=> no tageted reads statistic created)
  -T <file> Name for HTML Table row summary file (in output directory). Default: '' (=> none created)
  -l Log progress to STDERR. (A few primary progress messages will always be output.)
  -x Do not create the HTML file linking to all results created."

# should scan all args first for --X options
if [ "$1" = "--help" ]; then
    echo -e "$DESCR\n$USAGE\n$OPTIONS" >&2
    exit 0
fi

SHOWLOG=0
BEDFILE=""
GENOME=""
WORKDIR="."
OUTFILE=""
MAKEHML=1
RESHTML=""
ROWHTML=""
PADBED=""
BLOCKFILE=""
NONDUPREADS=0
UNIQUEREADS=0
ANNOBED=""
AMPOPT=""
TRIMP=0
PADVAL=0
TRGSID=""
RNABED=0
TRACKINGBED=""
FILTBAM=""
ISOLATE=0

while getopts "hladrtuxip:A:B:C:M:G:D:F:X:O:R:S:T:P:Q:" opt
do
  case $opt in
    A) ANNOBED=$OPTARG;;
    B) BEDFILE=$OPTARG;;
    C) TRGSID=$OPTARG;;
    D) WORKDIR=$OPTARG;;
    F) FILTBAM=$OPTARG;;
    G) GENOME=$OPTARG;;
    O) OUTFILE=$OPTARG;;
    P) PADBED=$OPTARG;;
    Q) BLOCKFILE=$OPTARG;;
    R) RESHTML=$OPTARG;;
    S) TRACKINGBED=$OPTARG;;
    T) ROWHTML=$OPTARG;;
    p) PADVAL=$OPTARG;;
    a) AMPOPT="-a";;
    d) NONDUPREADS=1;;
    r) RNABED=1;;
    t) TRIMP=1;;
    u) UNIQUEREADS=1;;
    x) MAKEHML=0;;
    i) ISOLATE=1;;
    l) SHOWLOG=1;;
    h) echo -e "$DESCR\n$USAGE\n$OPTIONS" >&2
       exit 0;;
    \?) echo $USAGE >&2
        exit 1;;
  esac
done
shift `expr $OPTIND - 1`

if [ $# -ne 2 ]; then
  echo "$CMD: Invalid number of arguments." >&2
  echo -e "$USAGE\n$OPTIONS" >&2
  exit 1;
fi

REFERENCE=$1
BAMFILE=$2

if [ -z "$GENOME" ]; then
  GENOME="$REFERENCE.fai"
fi
if [ -z "$RESHTML" ]; then
  RESHTML="results.html"
fi

BASECOVERAGE=1
if [ $RNABED -eq 1 ]; then
  AMPOPT="-r"
  BASECOVERAGE=0
fi

LOGOPT=''
if [ $SHOWLOG -eq 1 ]; then
  echo -e "\n$CMD BEGIN:" `date` >&2
  LOGOPT='-l'
fi
RUNPTH=`readlink -n -f $0`
RUNDIR=`dirname $RUNPTH`
if [ $SHOWLOG -eq 1 ]; then
  echo -e "RUNDIR=$RUNDIR\n" >&2
fi

if [ "$FILTBAM" = '-' ];then
  FILTBAM=''
fi

# Check environment

if ! [ -d "$RUNDIR" ]; then
  echo "ERROR: Executables directory does not exist at $RUNDIR" >&2
  exit 1;
elif ! [ -d "$WORKDIR" ]; then
  echo "ERROR: Output work directory does not exist at $WORKDIR" >&2
  exit 1;
elif ! [ -f "$GENOME" ]; then
  echo "ERROR: Genome (.fai) file does not exist at $GENOME" >&2
  exit 1;
elif ! [ -f "$REFERENCE" ]; then
  echo "ERROR: Reference sequence (fasta) file does not exist at $REFERENCE" >&2
  exit 1;
elif ! [ -f "$BAMFILE" ]; then
  echo "ERROR: Mapped reads (bam) file does not exist at $BAMFILE" >&2
  exit 1;
elif [ -n "$BEDFILE" -a ! -f "$BEDFILE" ]; then
  echo "ERROR: Reference targets (bed) file does not exist at $BEDFILE" >&2
  exit 1;
elif [ -n "$PADBED" -a ! -f "$PADBED" ]; then
  echo "ERROR: Padded reference targets (bed) file does not exist at $PADBED" >&2
  exit 1;
fi

# Get absolute file paths to avoid link issues in HTML
WORKDIR=`readlink -n -f "$WORKDIR"`
REFERENCE=`readlink -n -f "$REFERENCE"`
#BAMFILE=`readlink -n -f "$BAMFILE"`
GENOME=`readlink -n -f "$GENOME"`

BAMROOT=`echo $BAMFILE | sed -e 's/^.*\///'`
BAMNAME=`echo $BAMROOT | sed -e 's/\.[^.]*$//'`
ROOTNAME="$WORKDIR/$BAMNAME"

FILTBAM=`echo $FILTBAM | sed -e 's/\.bam$//'`

# Short descript of read filters
RTITLE=""
if [ $TRIMP -eq 1 ]; then
  RTITLE="Trimmed"
fi
if [ $UNIQUEREADS -eq 1 ]; then
  RTITLE="${RTITLE} Uniquely Mapped"
fi
if [ $NONDUPREADS -eq 1 ]; then
  RTITLE="${RTITLE} Non-duplicate"
fi
if [ -z "$RTITLE" ]; then
  RTITLE="All Mapped"
fi
RTITLE="$RTITLE Reads"

if [ "$OUTFILE" == "-" ]; then
  OUTFILE=""
fi
if [ -z "$OUTFILE" ]; then
  OUTFILE="${BAMNAME}.stats.txt"
fi
LOUTFILE="$OUTFILE"
OUTFILE="$WORKDIR/$OUTFILE"

# Get other derived options
FILTOPTS=""
SAMVIEWOPT="-F 4"
if [ $NONDUPREADS -eq 1 ];then
  FILTOPTS="-d"
  SAMVIEWOPT="-F 0x404"
fi
if [ $UNIQUEREADS -eq 1 ];then
  FILTOPTS="$FILTOPTS -u"
  SAMVIEWOPT="$SAMVIEWOPT -q 1"
fi
# BEDFILE and ANNOFILE become the same if only one is defined
BEDOPT=''
if [ -n "$BEDFILE" ]; then
  BEDOPT="-B \"$BEDFILE\""
fi
if [ -n "$ANNOBED" ]; then
  ANNOBEDOPT="-B \"$ANNOBED\""
  if [ -z "$BEDFILE" ]; then
    BEDFILE=$ANNOBED
    BEDOPT=$ANNOBEDOPT
  fi
else
  # called scripts fail w/o specific annotation format
  ANNOBEDOPT=''
fi

#--------- End command arg parsing ---------

if [ $SHOWLOG -eq 1 ]; then
  echo "" >&2
fi

############ Pre-run set up stuf goes here #################

if [ $TRIMP -eq 1 ]; then
  echo "(`date`) Trimming reads to targets..." >&2
  BAMTRIM="${WORKDIR}/${BAMNAME}.trim.bam"
  RT=0
  PTRIMCMD="java -Xmx8G -cp ${DIRNAME}/TRIMP_lib -jar ${DIRNAME}/TRIMP.jar \"$BAMFILE\" \"$BAMTRIM\" \"$REFERENCE\" \"$ANNOBED\""
  if [ $SHOWLOG -gt 0 ]; then
    echo "\$ $PTRIMCMD" >&2
  fi
  eval "$PTRIMCMD" || RT=$?
  if [ $RT -ne 0 ]; then
    echo "WARNING: TRIMP failed..." >&2
    echo "\$ $PTRIMCMD" >&2
  elif [ $SHOWLOG -gt 0 ]; then
    echo "> ${BAMTRIM}" >&2
  fi
  if [ -e "$BAMTRIM" ]; then
    SINDX="samtools index \"$BAMTRIM\"" >&2
    eval "$SINDX" || RT=$?
    if [ $RT -ne 0 ]; then
      echo "WARNING: samtools index failed... Proceeding with pre-trimmed BAM file." >&2
      echo "\$ $SINDX" >&2
    else
      BAMFILE="$BAMTRIM"
      BAMNAME="${BAMNAME}.trim"
      BAMROOT="${BAMNAME}.bam"
      if [ $SHOWLOG -gt 0 ]; then
        echo "> ${BAMTRIM}.bai" >&2
      fi
    fi
  else
    echo "WARNING: No trimmed BAM file found. Proceeding with pre-trimmed BAM file." >&2
  fi
  echo "" >&2
fi

########### Capture Title & User Options to Summary Text #########

echo -e "AmpliCat Report\n" > "$OUTFILE"
REFNAME=`echo "$REFERENCE" | sed -e 's/^.*\///' | sed -e 's/\.[^.]*$//'`
echo "Reference (File): $REFNAME" >> "$OUTFILE"
if [ -n "$BEDFILE" ]; then
  if [ -n "$TRGSID" ]; then
    TRGNAME=$TRGSID
  else
    TRGNAME=`echo "$BEDFILE" | sed -e 's/^.*\///' | sed -e 's/\.[^.]*$//'`
  fi
  echo "Targeted Regions: $TRGNAME" >> "$OUTFILE"
  if [ $PADVAL -gt 0 ]; then
    echo "Target Padding: $PADVAL" >> "$OUTFILE"
  fi
fi
ALMNAME=`echo "$BAMNAME" | sed -e 's/\.trim$//'`
if [ -n "$TRACKINGBED" ]; then
  RTITLE="$RTITLE and SampleID Tracking"
fi
echo "Alignments: $ALMNAME" >> "$OUTFILE"
#echo -e "Using: $RTITLE\n" >> "$OUTFILE"

########### Reduce BED to just non-overlaping targets #########

if [ "$ISOLATE" -ne 0 ];then
  TARGETBED="$WORKDIR/nonovlp_targets.bed"
  NONOVLPS=`$RUNDIR/isolate_amplicons.pl -n "$ANNOBED" "$TARGETBED"`
  if [ -z "$NONOVLPS" ];then
    echo -e "\nFailed to run isolate_amplicons.pl"
    echo "\$ $RUNDIR/isolate_amplicons.pl -n \"$ANNOBED\" \"$TARGETBED\"" >&2
    exit 1
  fi
  echo "Isolated targets analyzed:  $NONOVLPS" >> "$OUTFILE"
else
  TARGETBED="$ANNOBED"
  NONOVLPS=`awk '$1 !~ /^track/ {++c} END {print c}' "$ANNOBED"`
  echo "Targets analyzed:  $NONOVLPS" >> "$OUTFILE"
fi

# Add some basic alignment statistics

MAPPED_READS=`samtools view -c $SAMVIEWOPT "$BAMFILE"`
TREADS=$MAPPED_READS
#echo "Number of mapped reads:         $MAPPED_READS" >> "$OUTFILE"
if [ -n "$TRACKINGBED" ]; then
  TRACKING_READS=`samtools view -c $SAMVIEWOPT -L "$TRACKINGBED" "$BAMFILE"`
  TRACKING_READS=`echo "$TRACKING_READS $MAPPED_READS" | awk '{if($2<=0){$1=0;$2=1}printf "%.2f", 100*$1/$2}'`
  echo "Percent sample tracking reads:  $TRACKING_READS%" >> "$OUTFILE"
fi
#if [ -n "$BEDOPT" ]; then
#  TREADS=`samtools view -c $SAMVIEWOPT -L "$TARGETBED" "$BAMFILE"`
#  FREADS=`echo "$TREADS $MAPPED_READS" | awk '{if($2<=0){$1=0;$2=1}printf "%.2f", 100*$1/$2}'`
#  echo "Number of reads on target:      $TREADS" >> "$OUTFILE"
#  echo "Percent reads on target:        $FREADS%" >> "$OUTFILE"
#  if [ -n "$TRACKINGBED" ]; then
#    echo "Percent sample tracking reads:  $TRACKING_READS%" >> "$OUTFILE"
#  fi
#  if [ -n "$PADBED" ]; then
#    TREADS=`samtools view -c $SAMVIEWOPT -L "$PADBED" "$BAMFILE"`
#    FREADS=`echo "$TREADS $MAPPED_READS" | awk '{if($2<=0){$1=0;$2=1}printf "%.2f", 100*$1/$2}'`
#    echo "Percent reads on padded target: $FREADS%" >> "$OUTFILE"
#  fi
#elif [ -n "$TRACKINGBED" ]; then
#  echo "Percent sample tracking reads:  $TRACKING_READS%" >> "$OUTFILE"
#fi
#echo "" >> "$OUTFILE"

########### Primary analysis scripts go here #############

SAMOPT=''
if [ -n "$FILTBAM" ];then
  SAMOPT="-F \"$FILTBAM.sam\""
fi

TARGETCOVFILE="$ROOTNAME.amplicat.xls"

AMPCAT="$RUNDIR/classify_amp_reads.pl $LOGOPT $SAMOPT -a -b \"$BAMFILE\" \"$TARGETBED\" > \"$TARGETCOVFILE\""
eval "$AMPCAT" >&2
if [ $? -ne 0 ]; then
  echo -e "\nFailed to run classify_amp_reads.sh."
  echo "\$ $AMPCAT" >&2
  exit 1
fi

############ Convert filtered SAM to BAM ############

if [ -n "$FILTBAM" ];then
  echo "(`date`) Converting filtered SAM to BAM..." >&2
  eval "samtools view -h -S -b \"$FILTBAM.sam\" > \"$FILTBAM.bam\" 2> /dev/null";
  if [ $? -ne 0 ]; then
    echo " - Failed to convert $FILTBAM.sam to BAM file."
    rm -f "$FILTBAM.bam"
  else
    rm -f "$FILTBAM.sam"
    echo "> $FILTBAM.bam" >&2
    eval "samtools sort \"$FILTBAM.bam\" \"$FILTBAM.sort\" 2> /dev/null"
    eval "mv \"$FILTBAM.sort.bam\" \"$FILTBAM.bam\" 2> /dev/null"
    eval "samtools index \"$FILTBAM.bam\" 2> /dev/null"
    if [ $? -ne 0 ]; then
      echo " - Failed to index $FILTBAM.bam." >&2
    else
      echo "> $FILTBAM.bam.bai" >&2
    fi
  fi
fi

############ Secondary analysis scripts go here ##########

AMPCAT="$RUNDIR/summerize_amp_reads.pl \"$TARGETCOVFILE\" >> \"$OUTFILE\""
eval "$AMPCAT" >&2
if [ $? -ne 0 ]; then
  echo -e "\nFailed to run summerize_amp_reads.sh."
  echo "\$ $AMPCAT" >&2
fi

############ Detailed report workup goes here ############

if [ $MAKEHML -eq 1 ]; then
  if [ $SHOWLOG -eq 1 ]; then
    echo "" >&2
  fi
  echo -e "(`date`) Creating HTML report..." >&2
  if [ -n "$ROWHTML" ]; then
    ROWHTML="-T \"$ROWHTML\""
  fi
  GENOPT="-g"
  if [ -n "$BEDFILE" ]; then
    GENOPT=""
  fi
  SIDOPT=""
  if [ -n "$TRACKINGBED" ]; then
    SIDOPT="-i"
  fi
  COVERAGE_HTML="COVERAGE_html"
  PTITLE=`echo $BAMNAME | sed -e 's/\.trim$//'`
  HMLCMD="$RUNDIR/amplicat_report.pl -s \"$RTITLE\" $AMPOPT $ROWHTML $GENOPT $SIDOPT -N \"$BAMNAME\" -t \"$PTITLE\" -D \"$WORKDIR\" \"$COVERAGE_HTML\" \"$LOUTFILE\""
  eval "$HMLCMD" >&2
  if [ $? -ne 0 ]; then
    echo -e "\nERROR: amplicat_report.pl failed." >&2
    echo "\$ $HMLCMD" >&2
  elif [ $SHOWLOG -eq 1 ]; then
    echo "> ${WORKDIR}/$COVERAGE_HTML" >&2
  fi

  # Block Summary
#  if [ -n "$BLOCKFILE" ]; then
#    HMLCMD="perl $RUNDIR/amplicat_block.pl -s \"$RTITLE\" $AMPOPT $GENOPT $SIDOPT -O \"$BLOCKFILE\" -D \"$WORKDIR\" -S \"$LOUTFILE\" \"$BAMFILE\""
#    eval "$HMLCMD" >&2
#    if [ $? -ne 0 ]; then
#      echo -e "\nERROR: amplicat_block.pl failed." >&2
#      echo "\$ $HMLCMD" >&2
#    elif [ $SHOWLOG -eq 1 ]; then
#      echo "> ${WORKDIR}/${BLOCKFILE}" >&2
#    fi
#  fi
#  if [ $SHOWLOG -eq 1 ]; then
#    echo "HTML report complete: " `date` >&2
#  fi

fi

# Create local igv session file
if [ -n "$TSP_LIBRARY" ]; then
  TRACKOPT=''
  if [ -n "$ANNOBEDOPT" ]; then
    ANNOBED=`echo $ANNOBED | sed -e 's/^.*\///'`
    TRACKOPT="-a \"$ANNOBED\""
  fi
  BAMFILE=`echo $BAMFILE | sed -e 's/^.*\///'`
  COVCMD="$RUNDIR/create_igv_link.py -r ${WORKDIR} -b ${BAMFILE} $TRACKOPT -g ${TSP_LIBRARY} -s igv_session.xml"
  eval "$COVCMD" >&2
  if [ $? -ne 0 ]; then
    echo -e "\nWARNING: create_igv_link.py failed." >&2
    echo "\$ $COVCMD" >&2
  elif [ $SHOWLOG -eq 1 ]; then
    echo "> igv_session.xml" >&2
  fi
fi

############

if [ $SHOWLOG -eq 1 ]; then
  echo -e "\n$CMD END:" `date` >&2
fi

