#!/usr/bin/perl

#--------- Begin command arg parsing ---------

(my $CMD = $0) =~ s{^(.*/)+}{};
my $DESCR = "Analyze a given set of AmpliSeq aligned reads file for various phenomina:
strand bias, non-unique reads, fwd/rev read coverage, target false priming, outer false priming,
stacked early reads ends and primer-dimers. Optionally, simulated false priming may be performed
if the targets file includes primer sequences (with uracil bases for cleavage).
Outer false priming is indicated if the number of reads overlapping the amplicons (ovlp_reads)
does not equal the number of reads assignable to the amplicon (amp_reads), i.e. within the target
ends extended by the -S option. All other output metrics are based on percentage reads counts
relative to assignable ampicon reads (amp_reads). The false priming percentage FP_pc is based on the
numbers of reads starting late and reads stopping early that had an adapter detected (ZA tag).
Other reads ending early contribute to the EE_pc statistic only if there are at least (10) reads at
the same position. The PD_pc statistic is based on the number of reads within (12) bases of one
end of the insert. This might include reads also counting towards FP_pc or EE_pc. If the predicted
false priming is assessed, the extra fields FP_exp_pc and EE_exp_pc are added as the percentages of
read positons for FP_pc and EE_pc reads matching predicted false priming with uracil cutting. 
See options below for more information. Output is to STDOUT.";
my $USAGE = "Usage:\n\t$CMD [options] <BAM file> <targets file>";
my $OPTIONS = "Options:
  -h ? --help Display Help information
  -a Assign reads assuming Amplicon/AmpliSeq targets. (See -AD/-AU options below.)
     Default: Assign reads by maximum overlap with target. (Max. overlap also used to assign -a ties.)
  -b Targets file is in BED format, with no primer sequences (only 1st 4 BED columns read; -o option assumed).
     Default: File is tsv file with header and fields: chrom,start(1-base),end,amp_id,fprm_Useq,rprm_Useq,insert_seq.
  -d Ignore (PCR) Duplicate reads.
  -u Include only Uniquely mapped reads (MAPQ > 1).
  -f Use Full amplicon (with primers) as targets for false priming. Default: Use just insert regions.
  -l Print Log output to STDERR (including % progress and end-of-run statistics).
  -o Output only the (assumed) observed evidence for amplicon false priming. (-L,-P,-T and -U options ignorred.)
  -p Output only the predicted primer hybridizations for amplicon false priming. (Overrides -o and -P options.)
     Each primer hybridizing to each amplicon is output as separate record (line). <BAM file> argument may be ommited.
  -s Do not include reads that are soft-clipped (at 3') for early end counting.
  -D <N> Count reads that end only this many bases in to the target as possible primer-Dimers. Default: 20.
  -F <file> Create a Filtered SAM file for on-target reads not classified as produced by false priming. Default: ''.
  -L <N> Minimum 3' hybridization Length for false priming consideration. Minimum: 4. Default: 10.
  -O <N> Offset for read ends vs. insert ends to be consided possibly misprimed or short. Default: 5.
  -P <N> Add this number of columns of best primer predictions. Default: 0.
  -S <N> Ignore reads that Start more than this many bases before the target. Default: 30.
  -T <N> Minimum overlap Tm (4+2) cutoff considered for 3' hybridization. Default: 0.
  -U <N> Number of Uracil positions to cleave at from 3' (0,1,2). Default: 2.
  -AD <int> Amplicon Downstream limit for matching read start to target end (appropriate for +/- strand mapping).
     This assignment parameter is only employed if the -a option is provided. Default: 5.
  -AU <int> Amplicon Upstream limit for matching read start to target end (appropriate for +/- strand mapping).
     This assignment parameter is only employed if the -a option is provided. Default: 30.
  -AN <int> (Algorithm) Minimum Number of merged targets to use per samtools command. Default: 50.
     Forces more regions to be grouped than pysically overlapped within the -O option distance. The combination
     of -N, -O and -S option values affect run-time performance, depending mainly on the distribution of targets.
  -AO <int> (Algorithm) Minimum merged region separation (Overhang). Default 10000.
     The minimum base distance between one merged group and the next. Because reads can overlap the ends of
     multiple (merged) target regions, this value is forced to be at least 1000 to prevent a read being counted twice.
  -AS <int> (Algorithm) Maximum merged region Separation. Default 1000000. (Ineffective if -N option is < 2.)
     Limits the spacing between grouped merged target regions to reduce number of off-target reads sampled.
     This value is forced to be at least 10000.";

my $bedTargets = 0;
my $logopt = 0;
my $useSoftClipReads = 1;
my $fpSrtOffset = 5;
my $pdEndOffset = 20;
my $saSrtOffset = 30;

my $noFpp = 0;
my $outputFppOnly = 0;
my $fppUseFullAmp = 0;
my $fppTmCutoff = 0;
my $fppLenCutoff = 10;
my $fppNumUcuts = 2;
my $fppIgnoreUpos = 0;
my $fpp2ndUpos = 0;
my $fppNumSeqOutput = 0;

my $nondupreads = 0;
my $uniquereads = 0;

my $ampAlign = 0;
my $dsLimit = 5;
my $usLimit = 30;

my $minNumMerge = 50;
my $endOvlp = 10000;
my $maxMrgSep = 1000000;

my $samout = '';

# record stats for predicted FPs matching this fraction of observed FPs
my $highMatchFrac = 0.8; 
my $highMatchCount = 100; 

my $help = (scalar(@ARGV) == 0);
while( scalar(@ARGV) > 0 )
{
  last if($ARGV[0] !~ /^-/);
  my $opt = shift;
  if   ($opt eq '-a') {$ampAlign = 1;}
  elsif($opt eq '-b') {$bedTargets = 1;}
  elsif($opt eq '-d') {$nondupreads = 1;}
  elsif($opt eq '-u') {$uniquereads = 1;}
  elsif($opt eq '-f') {$fppUseFullAmp = 1;}
  elsif($opt eq '-l') {$logopt = 1;}
  elsif($opt eq '-o') {$noFpp = 1;}
  elsif($opt eq '-s') {$useSoftClipReads = 0;}
  elsif($opt eq '-p') {$outputFppOnly = 1;}
  elsif($opt eq '-D') {$pdEndOffset = int(shift);}
  elsif($opt eq '-F') {$samout = shift;}
  elsif($opt eq '-L') {$fppLenCutoff = int(shift);}
  elsif($opt eq '-O') {$fpSrtOffset = int(shift);}
  elsif($opt eq '-P') {$fppNumSeqOutput = int(shift);}
  elsif($opt eq '-S') {$saSrtOffset = int(shift);}
  elsif($opt eq '-T') {$fppTmCutoff = int(shift);}
  elsif($opt eq '-U') {$fppNumUcuts = int(shift);}
  elsif($opt eq '-AD') {$dsLimit = int(shift);}
  elsif($opt eq '-AU') {$dsLimit = int(shift);}
  elsif($opt eq '-AN') {$minNumMerge = int(shift);}
  elsif($opt eq '-AO') {$endOvlp = int(shift);}
  elsif($opt eq '-AS') {$maxMrgSep = int(shift);}
  elsif($opt eq '-h' || $opt eq "?" || $opt eq '--help') {$help = 1;}
  else
  {
    print STDERR "$CMD: Invalid option argument: $opt\n";
    print STDERR "$OPTIONS\n";
    exit 1;
  }
}
my $nargs = scalar(@ARGV);
my $noBam = ($nargs == 1 && $outputFppOnly);
if( $help )
{
  print STDERR "$DESCR\n";
  print STDERR "$USAGE\n";
  print STDERR "$OPTIONS\n";
  exit 1;
}
elsif( $nargs != 2 && !$noBam )
{
  print STDERR "$CMD: Invalid number of arguments.";
  print STDERR "$USAGE\n";
  exit 1;
}

my $bamfile = $noBam ? "" : shift;
my $ampfile = shift;

$fppLenCutoff = 4 if( $fppLenCutoff < 4 );
$fpp2ndUpos = 0 if( $fppIgnoreUpos );

# option dependencies and overrides
if( $bedTargets ) {
  $outputFppOnly = 0;
  $noFpp = 1;
}
my $outputFPexp = $outputFppOnly;
$noFpp = 0 if( $outputFppOnly );
$fppNumSeqOutput = 0 if( $noFpp );
$fppIgnoreUpos = 1 if( $fppNumUcuts <= 0 );
$fpp2ndUpos = 1 if( $fppNumUcuts >= 2 );

my $pdMaxLen = $saSrtOffset + $pdEndOffset;

$samout = '' if( $samout eq '-' );
my $samfilt = ($samout ne '' && $noBam == 0);

$minNumMerge = 1 if( $minNumMerge < 1 );
$endOvlp = 1000 if( $endOvlp < 1000 );
$maxMrgSep = 1000000000 if( $maxMrgSep < 10000 );
$maxMrgSep -= $endOvlp;

$usLimit *= -1;  # more convenient for testing

my $samopt= ($nondupreads ? "-F 0x704" : "-F 0x304").($uniquereads ? " -q 1" : "");

#--------- End command arg parsing ---------

# Read assays: amplicons and primers - arrays and $num_amps global to subroutines

open(AMPFILE,$ampfile) || die "Cannot open primers file $ampfile";
my $num_amps = $bedTargets ? 0 : -1;
my ($chrid,$srt,$end,$ampid,$fprimer,$rprimer,$aseq) = ("",0,0,"","","","");
while(<AMPFILE>) {
  # skip header
  if( $bedTargets ) {
    next if( /^track / );
    ++$num_amps;
  } else {
    next if( ++$num_amps == 0 );
  }
  chomp;
  if( $bedTargets ) {
    ($chrid,$srt,$end,$ampid) = split('\t',$_);
    $ampsrt[$num_amps] = $srt+1;
  } else {
    ($chrid,$srt,$end,$ampid,$fprimer,$rprimer,$aseq) = split('\t',$_);
    $ampsrt[$num_amps] = $srt+0;
  }
  $ampchr[$num_amps] = $chrid;
  $ampend[$num_amps] = $end+0;
  # pre-process the primers for Uracil positions
  $fupos[$num_amps] = $fppIgnoreUpos ? 0 : uCutPosition($fprimer);
  $rupos[$num_amps] = $fppIgnoreUpos ? 0 : uCutPosition($rprimer);
  if( $fpp2ndUpos ) {
    # these are offsets to the first U location
    $fupos2[$num_amps] = uCutPosition(substr($fprimer,0,$fupos[$num_amps])) - $fupos[$num_amps];
    $rupos2[$num_amps] = uCutPosition(substr($rprimer,0,$rupos[$num_amps])) - $rupos[$num_amps];
  }
  $ampid[$num_amps] = $ampid;
  $fprimer =~ tr/U/T/;
  $rprimer =~ tr/U/T/;
  $fprm[$num_amps] = $fprimer;
  $rprm[$num_amps] = $rprimer;
  $ampseq[$num_amps] = $fppUseFullAmp ? $fprimer.$aseq.revcomp($rprimer) : $aseq;
}
close(AMPFILE);
print STDERR "Loaded $num_amps amplicons\n" if( $logopt );

# Amplicon Group storage variables & arrays - global to subroutines
my ($ampGroupSrt,$ampGroupEnd,$lastNampRead) = (0,0,0);
my (@targOvpReads,@targReadFlgs,@targReadSrts,@targReadEnds,@targReadScrs,@targSamReads);
my (@ampReadData,$ampSamRead,$ampTotalOvlps);

# print header line, depending output options
my $xColPred = $fpp2ndUpos ? "" : "\n";
my $xColPFnd = ($fppNumSeqOutput > 0) ? "" : "\n";
if( $outputFppOnly ) {
  printf "amp_id\tlocus\tfp_pos\tfp_cut\tprm_id\tprm_len\tprm_tm$xColPred";
  print "\tfp_pos2\n" if( $fpp2ndUpos );
} else {
  my $extra = $noFpp ? "\n" : "\tFP_exp_pc\tEE_exp_pc$xColPFnd";
  print "amp_id\tlocus\tovlp_reads\tasgn_reads\tamp_reads\tfwd_pc\tfwd_cov_pc\trev_cov_pc\tNONU_pc\tMQ10_pc\tFP_pc\tEE_pc\tPD_pc\tsht_reads$extra";
  if( $fppNumSeqOutput == 1 ) {
    print "\tBestFPP\n";
  } elsif( $fppNumSeqOutput > 0 ) {
    for( my $i = 1; $i <= $fppNumSeqOutput; ++$i ) {print "\tBestFPP#$i"}
    print "\n";
  }
}

$|++; # turn on autoflush

# Progress updating and end-of-run statistics collection variables
my ($pcCmpNext,$pcCmpAdd) = (5,5);
my ($lowestTmHighmatch,$lowestTmFound,$lowestTmPossible,$highestTmNotFound) = (999,999,999,0);

# Hash to collect false priming matches at a position $fprecord{$pos} (per target amplicon)
# Values are array of (hyb_len, hyb_Tm, amp_idx, align_flgs). align_flgs:
#  b0 = targets fwd strand (3' cut), b1 = targets fwd strand (3' cut),
#  b2 = fwd/rev primer of amp idx,   b3 = 1st/2nd U position
my %fprecord;  # collects false priming matches at a position $fprecord{$pos}

# open output filtered SAM file and write out header
if( $samfilt ) {
  if( system("samtools view -H \"$bamfile\" > \"$samout\"") ) {
    print STDERR "Warning: Could not open bamfile to filter reads to $samout.\n";
    $samfilt = 0;
  } else {
    open( SAMOUT, ">>$samout" ) || die "Unexpected failure to open temporary SAM header file: $samout\n";
  }
}

# loop over amplicon regions and extract all read using samtools
for( my $ampn = 1; $ampn <= $num_amps; ++$ampn ) {
  # Analysze prediced FP positions
  if( !$noFpp ) {
    # Generates prediced only output, if flgged
    analyzeFalsePrime($ampn);
    if( $outputFppOnly ) {
      updateProgress($ampn) if( $logopt );
      next;
    }
  }
  # Read range limits relative to amplicon ends
  my $alen  = $ampend[$ampn] - $ampsrt[$ampn] + 1;
  my $asrt  = $ampsrt[$ampn] + $fpSrtOffset;
  my $aend  = $ampend[$ampn] - $fpSrtOffset;
  my $fpde  = $ampsrt[$ampn] + $pdEndOffset;
  my $rpde  = $ampend[$ampn] - $pdEndOffset;
  my $sasrt = $ampsrt[$ampn] - $saSrtOffset;
  my $saend = $ampend[$ampn] + $saSrtOffset;
  
  # Analyze reads for 'observed' FP locations
  my ($nreads,$freads,$rreads,$sreads,$covfwd,$covrev) = (0,0,0,0,0,0);
  my ($nmapq0,$nmapq10,$npdfwd,$npdrev) = (0,0,0,0);
  my (%foundFPs,%foundFPEnds,%foundPAs,%foundPAEnds);
  my ($lateSrts,$earlyEndsWithZATag,$earlyEndsNotLateSrts) = (0,0,0);
  my $samrec = 0;
  my @samrecs;
  while( getAmpRead($ampn) ) {
    my ($flg,$srt,$end,$scr,$zat) = @ampReadData;
    ++$nreads;
    $samrec = $nreads if( $samfilt );
    # skip analysis unless read is within limits of amplicon insert (=> not outer primed)
    # NOTE: outer primed reads are dis-counted in filtered reads as false primed
    next if( $srt < $sasrt || $end > $saend );
    my $keep_sclip = $useSoftClipReads || !($flg & 4);
    ++$sreads if( $keep_sclip && ($srt > $asrt || $end < $aend ));
    # only count early read ends if the read is not already a late start
    if( $flg & 16 ) {
      ++$rreads;
      if( $end <= $aend ) {
        ++$lateSrts;
        ++$foundFPs{$end};       # False prime by late start
        $foundFPEnds{$end} |= 2; # missing 5' sequence
        $samrec = 0;
      } elsif( $srt >= $asrt ) {
        # read start plus ZA length short vs. the fuzzy insert start => false primed
        if( $end-$zat > $asrt ) {
          ++$earlyEndsWithZATag;
          ++$foundFPs{$srt};       # False prime by ZA/early end
          $foundFPEnds{$srt} |= 1; # missing 3' sequence
          # test if this could also be classified as a primer-dimer
          ++$npdrev if( $end-$srt < $pdMaxLen );
          $samrec = 0;
        } elsif( $keep_sclip ) {
          ++$earlyEndsNotLateSrts;
          ++$foundPAs{$srt};       # trimmed sequence
          $foundPAEnds{$srt} |= 1; # missing 3' sequence
        }
      }
    } else {
      ++$freads;
      if( $srt >= $asrt ) {
        ++$lateSrts;
        ++$foundFPs{$srt};       # False prime by late start
        $foundFPEnds{$srt} |= 1; # missing 3' sequence
        $samrec = 0;
      } elsif( $end <= $aend ) {
        # read start plus ZA length short vs. the fuzzy insert end => false primed
        if( $srt+$zat <= $aend ) {
          ++$earlyEndsWithZATag;
          ++$foundFPs{$end};       # False prime by ZA/early end
          $foundFPEnds{$end} |= 2; # missing 5' sequence
          # test if this could also be classified as a primer-dimer
          ++$npdfwd if( $end-$srt < $pdMaxLen );
          $samrec = 0;
        } elsif( $keep_sclip ) {
          ++$earlyEndsNotLateSrts;
          ++$foundPAs{$end};       # trimmed sequence
          $foundPAEnds{$end} |= 2; # missing 5' sequence
        }
      }
    }
    # uniqueness tracking
    ++$nmapq0 if( $scr == 0 );
    ++$nmapq10 if( $scr <= 10 );
    # adjust ends to inserts for coverage calculations
    $srt = $ampsrt[$ampn] if( $srt < $ampsrt[$ampn] );
    $end = $ampend[$ampn] if( $end > $ampend[$ampn] );
    if( $flg & 16 ) {
      $covrev += $end - $srt + 1;
    } else {
      $covfwd += $end - $srt + 1;
    }
    # record possible reads for output (if not late starts or early ends for ZA tag reads)
    push( @samrecs, $samrec ) if( $samrec );
  }
  # Output sam reads for those that were not false priming
  printSamRecs( $ampn, \@samrecs ) if( $samfilt );

  # Calculate percentage stats from counts
  my $treads = $freads + $rreads;
  my $fpreads = $lateSrts + $earlyEndsWithZATag;
  my $pdreads = $npdfwd + $npdrev;
  my ($fwdpc,$fppc,$eepc,$pdpc) = (0,0,0,0);
  if( $treads > 0 ) {
    $fwdpc = 100 * $freads / $treads;
    $covfwd *= 100 / ($freads * $alen) if( $freads * $alen > 0 );
    $covrev *= 100 / ($rreads * $alen) if( $rreads * $alen > 0 );
    $nmapq0 *= 100 / $treads;
    $nmapq10 *= 100 / $treads;
    $fppc = 100 * $fpreads / $treads;
    $eepc = 100 * $earlyEndsNotLateSrts / $treads;
    $pdpc = 100 * $pdreads / $treads;
  }
  # Output for observed only output
  if( $noFpp ) {
    printf "%s\t%s:%d-%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n",
      $ampid[$ampn], $ampchr[$ampn], $ampsrt[$ampn], $ampend[$ampn], $ampTotalOvlps, $nreads, $treads,
      $fwdpc, $covfwd, $covrev, $nmapq0, $nmapq10, $fppc, $eepc, $pdpc, $sreads;
      updateProgress($ampn) if( $logopt );
      next;
  }
  
  # Determine number of reads explained by false priming predictions
  my ($tm,$ary);
  my ($fullExp,$possExp) = (0,0);
  my (%fppCnt,%fppDetail);
  while( ($pos,$cnt) = each(%foundFPs) ) {
    next if( !defined($fprecord{$pos}) );
    # ($hlen,$tm,$prmn,$rev)...
    my @fpr = @{$fprecord{$pos}};
    my $cntOnce = $cnt;
    for( my $i = 0; $i < scalar(@fpr); $i += 4 ) {
      # check for compatible end clipping
      next unless( defined($foundFPEnds{$pos}) && ($fpr[$i+3] & $foundFPEnds{$pos}) );
      $fullExp += $cntOnce;
      $cntOnce = 0;
      $prmid = ($fpr[$i+3] & 4 ? 'R-' : 'F-').$ampid[ $fpr[$i+2] ];
      $tm = $fpr[$i+1];
      $fppCnt{$prmid} += $cnt;
      push( @{$fppDetail{$prmid}}, ($cnt,$tm,$fpr[$i]) );
      $lowestTmFound = $tm if( $nreads >= $highMatchCount && $tm < $lowestTmFound );
    }
  }
  # Add in the possible FP matches
  while( ($pos,$cnt) = each(%foundPAs) ) {
    next if( !defined($fprecord{$pos}) );
    my @fpr = @{$fprecord{$pos}};
    my $cntOnce = $cnt;
    for( my $i = 0; $i < scalar(@fpr); $i += 4 ) {
      # check for compatible end clipping
      next unless( $fpr[$i+3] & $foundPAEnds{$pos} );
      $possExp += $cntOnce;
      $cntOnce = 0;
      $prmid = ($fpr[$i+3] & 4 ? 'R-' : 'F-').$ampid[ $fpr[$i+2] ];
      $tm = $fpr[$i+1];
      $fppCnt{$prmid} += $cnt;
      push( @{$fppDetail{$prmid}}, ($cnt,$tm,$fpr[$i]) );
      $lowestTmPossible = $tm if( $nreads >= $highMatchCount && $tm < $lowestTmPossible );
    }
  }
  # Add in remaining non-seen predicions, with fake hit count based on Tm (and by hyb. length)
  # Also resolve order for count matches by Tm for hits and determine lowest Tm for high match
  my $highMatchThres = $nreads >= $highMatchCount ? $highMatchFrac * ($fpreads + $earlyEndsNotLateSrts) : $highMatchCount;
  while( ($pos,$ary) = each(%fprecord) ) {
    my @fpr = @$ary;
    for( my $i = 0; $i < scalar(@fpr); $i += 4 ) {
      $prmid = ($fpr[$i+3] & 2 ? 'R-' : 'F-').$ampid[ $fpr[$i+2] ];
      $tm = $fpr[$i+1];
      if( defined($fppDetail{$prmid}) ) {
        $lowestTmHighmatch = $tm if( $fppCnt{$prmid} >= $highMatchThres && $tm < $lowestTmHighmatch );
        $fppCnt{$prmid} -= 1/$tm;
        next;
      }
      $cnt  = -1/$tm - 0.001/$fpr[$i]; # hyb. length should just split ties
      $fppCnt{$prmid} += $cnt;
      push( @{$fppDetail{$prmid}}, (0,$tm,$fpr[$i]) );
      $highestTmNotFound = $tm if( $nreads >= $highMatchCount && $tm > $highestTmNotFound );
    }
  }
  # Output for supporting predicted false priming...
  my $fp_exp_pc = $fpreads ? 100 * $fullExp / $fpreads : 0;
  my $ee_exp_pc = $earlyEndsNotLateSrts ? 100 * $possExp / $earlyEndsNotLateSrts : 0;
  printf "%s\t%s:%d-%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d$xColPFnd",
    $ampid[$ampn], $ampchr[$ampn], $ampsrt[$ampn], $ampend[$ampn], $ampTotalOvlps, $nreads, $treads,
    $fwdpc, $covfwd, $covrev, $nmapq0, $nmapq10, $fppc, $eepc, $pdpc, $fp_exp_pc, $ee_exp_pc, $sreads;

  # Add the top primer FP prediction hits to output
  if( $fppNumSeqOutput > 0 ) {
    my $numOut = 0;
    foreach my $prmid (sort {$fppCnt{$b} <=> $fppCnt{$a}} keys %fppCnt) {
      printf "\t%s (%s)", $prmid, join(',',@{$fppDetail{$prmid}});
      last if( ++$numOut >= $fppNumSeqOutput );
    }
    # fill in columns if not enough predicted fps
    print "\tNA" while( ++$numOut <= $fppNumSeqOutput );
    print "\n";
  }
  updateProgress($ampn) if( $logopt );
}
close(SAMOUT) if( $samfilt );

# end-of-run summary
if( $logopt ) {
  updateProgress($num_amps);
  if( !$noFpp ) {
    print STDERR "For amplicons having at least $highMatchCount reads:\n";
    $lowestTmHighmatch = '--' if( $lowestTmHighmatch == 999 );
    $lowestTmFound = '--' if( $lowestTmFound == 999 );
    $lowestTmPossible = '--' if( $lowestTmPossible == 999 );
    $highestTmNotFound = '--' if( $highestTmNotFound == 0 );
    printf STDERR " Lowest predicted FP Tm explaining %.0f%%+ observations: $lowestTmHighmatch\n",100*$highMatchFrac;
    print STDERR " Lowest predicted FP Tm matching strong observations: $lowestTmFound\n";
    print STDERR " Lowest predicted FP Tm matching weak observations:   $lowestTmPossible\n";
    print STDERR " Highest predicted FP Tm matching no observations:    $highestTmNotFound\n";
  }
}

# ----------------- End ------------------

sub updateProgress
{
  my $ampn = $_[0];
  my $pcCmp = 0;
  # delay printing upto the greatest amount completed - for small input data sets
  while( 100*$ampn/$num_amps >= $pcCmpNext ) {
    $pcCmp = $pcCmpNext;
    $pcCmpNext += $pcCmpAdd;
  }
  printf STDERR "%3d%% complete.\n",$pcCmp if( $pcCmp > 0 );
}

sub uCutPosition
{
  my $up = rindex($_[0],'U');
  return $up < 0 ? 0 : $up;
}

sub revcomp
{
  my $seq = $_[0];
  $seq =~ tr/ACGTacgt/TGCAtgca/;
  return scalar reverse $seq;
}

sub analyzeFalsePrime
{
  my $q = $_[0];
  my $ampfwd = $ampseq[$q];
  %fprecord = ();  # reset w/o freeing (all) memory
  my $amprev = revcomp($ampfwd);
  for( my $p = 1; $p <= $num_amps; ++$p ) {
    # only ignore self-priming if not looking at full amplicon -> fix later
    next if( $p == $q && $fppUseFullAmp );
    storeFPRecord($q,$p,$ampfwd,0);
    storeFPRecord($q,$p,$amprev,1);
  }
}

sub storeFPRecord
{
  # same as storeFPRecord but extra formatted printing
  # - not so easy to get from the extra hash elements created/added to!
  my ($ampn,$prmn,$amps,$rev) = @_;
  my $end = $rev ? 'R' : 'F';
  my $deb = $rev ? 2 : 1;
  my ($hlen,$hpos,$tm) = maxovlpAboveTm($fprm[$prmn],$amps);
  if( $hlen > 0 ) {
    # position is end of the read, i.e. one base to the left or right of the U cut position
    $hpos = $rev ? $ampend[$ampn]-$hpos-$fupos[$prmn]-1 : $ampsrt[$ampn]+$hpos+$fupos[$prmn]+1;
    push( @{$fprecord{$hpos}}, ($hlen,$tm,$prmn,$deb) );
    print "$ampid[$ampn]\t$ampchr[$ampn]:$ampsrt[$ampn]-$ampend[$ampn]\t$hpos\t$end\tF-$ampid[$prmn]\t$hlen\t$tm$xColPred" if( $outputFPexp );
    if( $fpp2ndUpos ) {
      $hpos += $rev ? -$fupos2[$prmn] : $fupos2[$prmn];
      push( @{$fprecord{$hpos}}, ($hlen,$tm,$prmn,$deb+8) );    # end + 2nd U
      printf "\t$hpos\n" if( $outputFPexp );
    }
  }
  ($hlen,$hpos,$tm) = maxovlpAboveTm($rprm[$prmn],$amps);
  if( $hlen > 0 ) {
    $hpos = $rev ? $ampend[$ampn]-$hpos-$rupos[$prmn]-1 : $ampsrt[$ampn]+$hpos+$rupos[$prmn]+1;
    push( @{$fprecord{$hpos}}, ($hlen,$tm,$prmn,$deb+4) );  # end + R primer
    print "$ampid[$ampn]\t$ampchr[$ampn]:$ampsrt[$ampn]-$ampend[$ampn]\t$hpos\t$end\tR-$ampid[$prmn]\t$hlen\t$tm$xColPred" if( $outputFPexp );
    if( $fpp2ndUpos ) {
      $hpos += $rev ? -$rupos2[$prmn] : $rupos2[$prmn];
      push( @{$fprecord{$hpos}}, ($hlen,$tm,$prmn,$deb+12) );  # end + R primer + 2nd U
      printf "\t$hpos\n" if( $outputFPexp );
    }
  }
}

sub maxovlpAboveTm
{
  my ($qry,$trg) = @_;
  my ($hlen,$hpos) = maxovlp($qry,$trg);
  return (0,0,0) if( $hlen <= 0 );
  my $tm = tempGC( substr($qry,-$hlen) );
  return (0,0,0) if( $tm < $fppTmCutoff );
  return ($hlen,$hpos,$tm);
}

sub maxovlp
{
  my ($qry,$trg) = @_;
  my $qlen = length($qry);
  my ($bidx,$blen) = (0,0);
  my $idx;
  # look for smallest length match first: if this fails no point looking further!
  for( my $len = $fppLenCutoff; $len <= $qlen; ++$len )
  {
    # seach from 3' since this is the more likely observed hybdirdization
    $idx = rindex($trg,substr($qry,-$len));
    last if($idx < 0);
    $bidx = $idx;
    $blen = $len;
  }
  # return match position for 5' end of aligned primer (not at match start)
  return ($blen,$bidx+$blen-$qlen);
}

sub tempGC
{
  my $seq = $_[0];
  my $tm = 4*($seq =~ s/[cgGC]//g) + 2*($seq =~ s/[atAT]//g);
}

# Essentially an iterator function that loads a group of assigned on-target read for a given amplicon index and returns
# next read to global data values and 0/1 depending on whether a call was made beyond the last read for the given amplicon.
# Data for the read just made is stored in @ampReadData (and $ampSamRead if the filtered sam output is used).
# But the total number of overlapping read for the amplicon is only available after all reads have been iterated (on return of 0).
# Assumes only sequential fwd access, although can be used for random access by externally resetting $ampGroupSrt = 0.
#
sub getAmpRead
{
  my $namp = $_[0];
  if( $ampGroupSrt == 0 || $namp >= $ampGroupEnd ) {
    loadAmpGroup($namp);
    $lastNampRead = -1;
  }
  my $ampt = $namp - $ampGroupSrt;
  if( ++$lastNampRead >= scalar(@{$targReadSrts[$ampt]}) ) {
    $lastNampRead = -1;
    $ampTotalOvlps = $targOvpReads[$ampt];
    return 0;
  }
  @ampReadData = ( $targReadFlgs[$ampt][$lastNampRead], $targReadSrts[$ampt][$lastNampRead],
    $targReadEnds[$ampt][$lastNampRead], $targReadScrs[$ampt][$lastNampRead], $targReadZATs[$ampt][$lastNampRead] );
  $ampSamRead = $targSamReads[$ampt][$lastNampRead] if( $samfilt );
  return 1;
}

# Cannot stream on-target reads for a specific target one at a time because targets may overlap.
# Ideally the merged subgroups (for each samtools invocation) could be processed separately for the sake of minimal
# memory usage, but this would mean making the code re-enterant for each merged subgroup of targets and quite tricky.

sub loadAmpGroup
{
  # starting at the given amplicon index, create list of amplicons that are grouped for single samtools view retrieval
  $ampGroupSrt = $_[0];
  my $anum = $ampGroupSrt;
  my $achr = $ampchr[$anum];
  my $grpSrt = $ampsrt[$anum];
  my $mrgEnd = $ampend[$anum] + $endOvlp;
  my ($padend,$sepend,$samrec);
  my $ntrg = 1;
  while( ++$anum < $num_amps && $ampchr[$anum] eq $achr ) {
    # force addition of next target that is within merge range of last
    $sepend = $ampsrt[$anum] - $mrgEnd;
    if( $sepend <= 0 ) {
      $padend = $ampend[$anum]+$endOvlp;
      $mrgEnd = $padend if( $padend > $mrgEnd );
      ++$ntrg;
      next;
    }
    # do not add any more (merged) targets if next starts too far from last or minimum group size is met
    last if( $sepend >= $maxMrgSep || $ntrg >= $minNumMerge );
    $mrgEnd = $ampend[$anum]+$endOvlp;
    ++$ntrg;
  }
  $ampGroupEnd = $anum;

  # initialize retrived data storage
  @targOvpReads = (0) x $ntrg;
  @targReadFlgs = ();
  @targReadSrts = ();
  @targReadEnds = ();
  @targReadScrs = ();
  @targReadZATs = ();
  @targSamReads = ();

  # record all reads overlapping grouped proximal targets and assign to specific targets
  my $firstRegion = $ampGroupSrt;
  my $grpSrt = $ampsrt[$ampGroupSrt];
  my $grpEnd = $ampend[$anum-1];
  open( MAPPINGS, "samtools view $samopt \"$bamfile\" $achr:$grpSrt-$grpEnd |" )
    || die "Failed to pipe reads from $bamfile for regions in $bedfile\n";
  while( <MAPPINGS> ) {
    $samrec = $_;
    my $zat = /\tZA:i:(\d+)/ ? $1 : 9999;  # if no ZA tag, set value to way high so extra undefined test not required
    my ($rid,$flg,$chr,$srt,$scr,$cig) = split('\t',$_);
    $srt += 0;
    my $end = $srt-1;
    # check and add on soft-clip to aligned read end
    if( $cig =~ s/(\d+)S// ) {
      $flg & 16 ? $srt -= $1 : $end += $1;
      $flg |= 4; # hijack the 0x4 bit for indicating soft-clipping (unmapped read flag no longer relevant)
    }
    # determine alignment end
    while( $cig =~ s/^(\d+)(.)// ) {
      # here we are including soft-clips as part of the match for better read-length/alignment measurement but excluding inserts
      $end += $1 if( index("MDX=",$2) >= 0 );
    }
    # determine best target/amplicon overlap
    my ($maxOvlp,$maxOvlp2,$trgIndx,$trgIndx2,$noAmpEndAlign) = (-1,-1,0,0,0);
    for( $anum = $firstRegion; $anum < $ampGroupEnd; ++$anum ) {
      # safe to stop early when read end is prior to start of target
      my $asrt = $ampsrt[$anum];
      last if( $end < $asrt );
      # no match if read start is after target end, but target ends can overlap therefore need to carry on searching
      my $aend = $ampend[$anum];
      if( $srt > $aend ) {
        # adjust start of list for further reads if no earlier target end found
        $firstRegion = $anum+1 if( $maxOvlp < 0 );
        next;
      }
      # record a hit for all targets overlaped by this read
      $ntrg = $anum - $ampGroupSrt;
      ++$targOvpReads[$ntrg];
      my $dsrt = $srt - $asrt;
      my $dend = $aend - $end;
      # test if this can be assigned to an amplicon
      if( $ampAlign ) {
        my $d3pend = ($flg & 16) ? $dend : $dsrt;
        $noAmpEndAlign = ($d3pend < $usLimit || $d3pend > $dsLimit);
      }
      # chop off target overhang from for overlap assessment
      $dsrt = 0 if( $dsrt < 0 );
      $dend = 0 if( $dend < 0 );
      $asrt = $aend - $asrt - $dsrt - $dend; # actually 1 less than overlap
      # record separate maxima for amplicon end matched reads and others
      # - in case of a tie, keep the most 3' match for backwards-compatibility to old 3.6 version
      if( $noAmpEndAlign ) {
        if( $asrt >= $maxOvlp2 ) {
          $maxOvlp2 = $asrt;
          $trgIndx2 = $ntrg;
        }
      } elsif( $asrt >= $maxOvlp ) {
        $maxOvlp = $asrt;
        $trgIndx = $ntrg;
      }
    }
    # only use best result matched to any amplicon if there was no best end amplicon alignment was found
    if( $maxOvlp < 0 ) {
      $maxOvlp = $maxOvlp2;
      $trgIndx = $trgIndx2;
    }
    # save this read for assigned target (some reads may be between targets)
    if( $maxOvlp >= 0 ) {
      # Note: flg, srt and scr could easily be re-parsed from samrec.
      # Here they are not for greater peformance, particularly if the filtered sam option is not required.
      push( @{$targReadFlgs[$trgIndx]}, $flg );
      push( @{$targReadSrts[$trgIndx]}, $srt );
      push( @{$targReadEnds[$trgIndx]}, $end );
      push( @{$targReadScrs[$trgIndx]}, $scr );
      push( @{$targReadZATs[$trgIndx]}, $zat );
      push( @{$targSamReads[$trgIndx]}, $samrec ) if( $samfilt );
      
    }
  }
  close( MAPPINGS );
}

sub printSamRecs
{
  my $ampn = $_[0] - $ampGroupSrt;
  my $recs = $_[1];
  for( my $i = 0; $i < scalar(@$recs); ++$i ) {
    print SAMOUT $targSamReads[$ampn][$recs->[$i]-1] if( $recs->[$i] );
  }
}

