#!/usr/bin/perl
# Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved

use strict;
use warnings;

#--------- Begin command arg parsing ---------

(my $CMD = $0) =~ s{^(.*/)+}{};
my $DESCR = "Create down-sampling matrix to estimated 1x representation rate and/or target down-samping percentage
given a matirix of target (amplicon) reads for each strand for each barcode.";
my $USAGE = "Usage:\n\t$CMD [options] <BarcodeMatrix file> <DownSampleMatrix file>\n";
my $OPTIONS = "Options:
  -h ? --help Display Help information
  -d Down-sample to target level only. Do not down-sample to calculated relatative representation ratios.
  -l Log output to STDERR.
  -m Assumed Merged barcodes for down-sampling but with fwd and rev reads downsampled separately.
     For deduping by representation rate (to minimum -T level) total reads per amplicon will be the same but the initial
     barcode representation (number of reads) is ignorred and may be reduced to 0.
     For downsampling to target level (-D, post deduping) this means individual barcode coverage may be reduced to 0.
  -r Print end-of-run down-sampled Representation statistics to STDOUT.
  -s Print end-of-run Summary statistics to STDOUT.
  -t Output statistics in HTML format.
  -D <N> Down-sample to this number of total reads as a percentage of all reads in the input coverage matrix. Default: 100%.
  -F <file> Output target representation Fold (duplication) matrix to the given file.
  -P <N> Number of decimal Places to (pre)round fold representation values. Minimum: 0. Default: 1.
  -R <N> Number of fold Represntation bins output in summary statistics. Minimum: 3. Default: 5.
  -T <N> Target coverage to set threshold for duplication rate to downsample over-represented samples.
     This is percentage of initial read coverage. It is only applied if sampling to 1x representation is not already
     above this limit. Hence, 0 => use whatever the coverage is left after sampling down to 1x. Default: 0.";

# Acdording to en.wikipedia.org/wiki/Quartile:
# Method 1 (TI-83 calculator boxplot and "1-Var Stats" functions)
# Method 2 - gives higher Q1+Q2 values and seems to match MS Excel calculations
# - Other methods not suported
my $quartile_method = 2;

my $logstats = 0;
my $logds = 0;
my $logopt = 0;
my $repbins = 5; # must be at least 3
my $htmlstats = 0;
my $dpround = 1; # decimal places to round off to, e.g. 0 => nearest integer, 1 => 1 decimal place, etc.
my $dsonly = 0;
my $mergebcs = 0;
my $targetDSFrac = 100;

my $dedupTarget = 0;

my $subbc = 0;	# number of subset barcodes to process 0 => all in range
my $bcsrt = 1;  # first barcode to start with
my $bcinc = 1;  # every Nth barcode to take

my $dupmatirx = '';

my $num_dstries = 4;

my $help = (scalar(@ARGV) == 0);
while( scalar(@ARGV) > 0 )
{
  last if($ARGV[0] !~ /^-/);
  my $opt = shift;
  if   ($opt eq '-d') {$dsonly = 1;}
  elsif($opt eq '-l') {$logopt = 1;}
  elsif($opt eq '-m') {$mergebcs = 1;}
  elsif($opt eq '-r') {$logds = 1;}
  elsif($opt eq '-s') {$logstats = 1;}
  elsif($opt eq '-t') {$htmlstats = 1;}
  elsif($opt eq '-D') {$targetDSFrac = 0+shift;}
  elsif($opt eq '-F') {$dupmatirx = shift;}
  elsif($opt eq '-R') {$repbins = int(shift);}
  elsif($opt eq '-P') {$dpround = 0+shift;}
  elsif($opt eq '-T') {$dedupTarget = 0+shift;}
  elsif($opt eq '-h' || $opt eq "?" || $opt eq '--help') {$help = 1;}
  else
  {
    print STDERR "$CMD: Invalid option argument: $opt\n";
    print STDERR "$OPTIONS\n";
    exit 1;
  }
}
my $nargs = scalar(@ARGV);
if( $help )
{
  print STDERR "$DESCR\n";
  print STDERR "$USAGE\n";
  print STDERR "$OPTIONS\n";
  exit 1;
}
elsif( $nargs != 2 )
{
  print STDERR "$CMD: Invalid number of arguments.";
  print STDERR "$USAGE\n";
  exit 1;
}

my $bcmatrix = shift;
my $dsmatrix = shift;

$dupmatirx = '' if( $dsmatrix eq '-' );

my $nlchr = $htmlstats ? "<br/>\n" : "\n";
my $iform = $dpround < 0 ? "%g" : "%.${dpround}f";

$repbins = 3 if( $repbins < 3 );
$subbc = 0 if( $subbc < 0 );
$bcsrt = 0 if( --$bcsrt < 0 );
$bcinc = 1 if( $bcinc < 1 );

$targetDSFrac = 100 if( $targetDSFrac <= 0 || $targetDSFrac > 100 );
$dedupTarget = 0 if( $dedupTarget < 0 );
if( $dedupTarget >= 100 ) {
  $dedupTarget = 100;
  $dsonly = 1;
}

my $dupMatrix = ($dupmatirx ne '');

#--------- End command arg parsing ---------

print  STDERR "  Assuming merged barcode coverage.\n" if( $mergebcs );
unless( $dsonly ) {
  if( $dedupTarget > 0 ) {
    print  STDERR "  De-duping to >= $dedupTarget% target coverage.\n";
  } else {
    print  STDERR "  De-duping to 1x target coverage.\n";
  }
}
printf STDERR "  Down-sampling to %g%% target coverage.\n",$targetDSFrac if( $targetDSFrac < 100 );

open( TABLEIN, "$bcmatrix" ) || die "Cannot open barcode/amplicon coverage matrix file $bcmatrix";
open( DSMATOUT, ">$dsmatrix" ) || die "Cannot write to barcode/amplicon downsample matrix file $dsmatrix";
if( $dupMatrix ) {
  open( DUPMATRIX, ">$dupmatirx" ) || die "Cannot write to barcode/amplicon duplication matrix file $dupmatirx";
}

# set random seed using the bcmatrix file name
my $sdstr = $bcmatrix;
$sdstr =~ s/.*\///;
setRandSeed($CMD.$sdstr);

my (@rows,@ampid,@locus);
my ($numbcs, $linenum) = (0,0);
while(<TABLEIN>) {
  chomp;
  my @fields = split('\t');
  unless($linenum++) {
    # determine initial number of barcodes from first line and validate sub-barcode selection
    my $amp_id = shift(@fields);
    my $loc_id = shift(@fields);
    unless( $numbcs ) {
      $numbcs = scalar(@fields);
      if( $bcsrt > 0 || $bcinc > 1 ) {
        my $nsam = ($numbcs - $bcsrt) / $bcinc;
        $nsam = int($nsam) + ($nsam > int($nsam)); # ceil($nsam)
        $subbc = $nsam if( $subbc < $nsam );
      }
      $subbc = 0 if( $subbc > $numbcs );
    }
    # print adjusted headers
    print DUPMATRIX "$amp_id\t$loc_id";
    print DSMATOUT "$amp_id\t$loc_id";
    for( my $bc = $bcsrt; $bc < $numbcs; $bc += $bcinc ) {
      print DUPMATRIX "\t$fields[$bc]";
      print DSMATOUT "\t$fields[$bc]";
    }
    print DUPMATRIX "\n";
    print DSMATOUT "\n";
    next;
  }
  push( @ampid, shift(@fields) );
  push( @locus, shift(@fields) );
  # restrict to just the barcodes of interest
  if( $subbc ) {
    my @data;
    for( my $bc = $bcsrt; $bc < $numbcs; $bc += $bcinc ) {
      push( @data, $fields[$bc] );
    }
    push( @rows, \@data );
  } else {
    push( @rows, \@fields );
  }
}
close(TABLEIN);

my $numamps = scalar(@rows);
if( $numamps <= 0 ) {
  close(DUPMATRIX);
  close(DSMATOUT);
  print STDERR "$CMD: No amplicon was read from $bcmatrix.";
  exit 1;
}
# reset to the number of barcodes actually read
$numbcs = scalar(@{$rows[0]});

# make copy of original array for down-sizing
my @dsrows = map {[@$_]} @rows;

# determine median/mean stats amplicon representation per barcode
my @bcmedian;
my (@bcLLT,@bcULT);
for( my $bc = 0; $bc < $numbcs; ++$bc ) {
  my @data;
  for( my $an = 0; $an < $numamps; ++$an ) {
    push( @data, $rows[$an][$bc] );
  }
  my ($Q1,$Q2,$Q3) = quartiles(\@data);
  my $IQR = $Q3 - $Q1;
  my $LLT = $Q1 - $IQR;
  my $ULT = $Q1 - 0.5*$IQR;
  # adjust the lower range to 1, since we should have minimally one read (in both directions)
  $LLT = 1 if( $LLT < 1 );
  $ULT = $LLT+1 if( $ULT <= $LLT );
  printf STDERR "> Median amplicon coverage for barcode %d = %g\n",$bc+1,$Q2 if( $logopt );
  printf STDERR ">   Q1 = %g, Q3 = %g, IQR = %g, LLT-ULT = %g-%g\n",$Q1,$Q3,$IQR,$LLT,$ULT if( $logopt );
  push( @bcmedian, $Q2 );
  push( @bcLLT, $LLT );
  push( @bcULT, $ULT );
#  my ($av,$sd) = avedev(\@data);
#  push( @bcLLT, $av-2*$sd );
#  printf STDERR "> Ave/SD amplicon coverage for barcode %d = %g, %g\n",$bc+1,$av,$sd if( $logopt );
}

# determine the median absolute deviation across barcodes (before hiding zeros)
my $bcmd_md = median(\@bcmedian);
if( $logstats ) {
  my @absdev;
  for( my $bc = 0; $bc < $numbcs; ++$bc ) {
    push( @absdev, abs($bcmedian[$bc] - $bcmd_md) );
  }
  my $bcmd_mad = median(\@absdev);
  print "${nlchr}Amplicon representation analysis:${nlchr}${nlchr}";
  printf "Median of median amplicon coverage over 2x%d barcodes: %.1f$nlchr",($numbcs/2),$bcmd_md; 
  printf "   MAD of median amplicon coverage over 2x%d barcodes: %.1f$nlchr",($numbcs/2),$bcmd_mad; 
}

# set the minimum down-size coverage threshold, i.e. cannot reduce the number of initial reads below this value
my $minDSCov = int(0.5+0.1*$bcmd_md);  # 10% of median value
# always leave at least one read strand per per amplicon per barcode, if not already 0
$minDSCov = 1 if( $minDSCov < 1 );

# remove any zero medians in case any barcode had all amplicons at 0 coverage
for( my $bc = 0; $bc < $numbcs; ++$bc ) {
  $bcmedian[$bc] = 1 if( $bcmedian[$bc] <= 0 );
}

# count lower outliers and upper outliers per amplicon
my ($ampLLO,$ampULO) = (0,0);

# normalize all reads to median values for barcode
for( my $an = 0; $an < $numamps; ++$an ) { 
  my $row = $rows[$an];
  my ($lowerOutliers,$upperOutliers) = (0,0);
  for( my $bc = 0; $bc < $numbcs; ++$bc ) {
    ++$lowerOutliers if( $row->[$bc] < $bcLLT[$bc] );
    ++$upperOutliers if( $row->[$bc] < $bcULT[$bc] );
    $row->[$bc] /= $bcmedian[$bc];
  }
  if( $lowerOutliers ) {
    ++$ampLLO if( $lowerOutliers == $numbcs );
    ++$ampULO if( $upperOutliers == $numbcs ); # only if at least one lower outlier
    printf "> Amplicon %s has %d Q1-IQR & %d Q1-IQR/2 outliers$nlchr", $ampid[$an], $lowerOutliers, $upperOutliers if( $logopt );
  }
}

# for each amplicon re-normalize by the median of the non-zero normalized barcode representations
for( my $an = 0; $an < $numamps; ++$an ) {
  my $row = $rows[$an];
  my @data;
  for( my $bc = 0; $bc < $numbcs; ++$bc ) {
    push( @data, $row->[$bc] ) if( $row->[$bc] > 0 );
  }
  # if any non-0 values look for lowest non outlier and re-normalize
  if( scalar(@data) ) {
    my $med = median(\@data);
    for( my $bc = 0; $bc < $numbcs; ++$bc ) {
      $row->[$bc] /= $med;
      $row->[$bc] = sprintf($iform,$row->[$bc]);
    }
  }
}

if( $dupMatrix ) {
  # output relative amplicon expression levels
  for( my $an = 0; $an < $numamps; ++$an ) {
    my $row = $rows[$an];
    printf DUPMATRIX "%s\t%s", $ampid[$an], $locus[$an];
    for( my $bc = 0; $bc < $numbcs; ++$bc ) {
      printf DUPMATRIX "\t%g", $row->[$bc];
    }
    print DUPMATRIX "\n";
  }
  close(DUPMATRIX);
}

# analyze degree of over performing assays vs. normalzed representation
if( $logstats ) {
  my @narcnts = (0) x $repbins;
  my @nar2sds = (0) x $repbins;
  for( my $an = 0; $an < $numamps; ++$an ) {
    my $row = $rows[$an];
    my ($av,$sd) = avedev($row);
    my $sd2cut = $av + 2*$sd;
    my @nrbins = (0) x $repbins;
    my @nsbins = (0) x $repbins;
    for( my $bc = 0; $bc < $numbcs; ++$bc ) {
      for( my $ct = 1; $ct < $repbins; ++$ct ) {
        #++$nrbins[$ct] if( $row->[$bc] >= (1<<($ct-1))-0.5 );
        if( $row->[$bc] >= (1<<($ct-1))-0.5 ) {
          ++$nrbins[$ct];
          ++$nsbins[$ct] if( $row->[$bc] > $sd2cut );
          #print STDERR "Amplicon $ampid[$an] has a ${ct}x hit at $row->[$bc] rep. (AV = $av, SD = $sd, AV+2SD = $sd2cut)\n" if( $ct > 1 ); 
        }
      }
    }
    # adjust for non-cumulative stats
    $nrbins[0] = $numbcs - $nrbins[1];
    $nrbins[1] = (($nrbins[1]-$nrbins[2]) == $numbcs);
    $nsbins[0] = $numbcs - $nsbins[1];
    $nsbins[1] = (($nsbins[1]-$nsbins[2]) == $numbcs);
    # record hits per amplicon
    for( my $ct = 0; $ct < $repbins; ++$ct ) {
      ++$narcnts[$ct] if( $nrbins[$ct] );
      ++$nar2sds[$ct] if( $nsbins[$ct] );
    }
  }
  # print stats
  printf "%6.2f%% amplicons were occasional drop-outs or consistent low performers.$nlchr", 100*$ampULO/$numamps;
  printf "%6.2f%% amplicons were consistent drop-outs or lowest performers.$nlchr", 100*$ampLLO/$numamps;
  #printf "%6.2f%% amplicons have some barcodes at < 1x representation.$nlchr", 100*$narcnts[0]/$numamps;
  printf "%6.2f%% amplicons have ALL barcodes at ~ 1x representation.$nlchr", 100*$narcnts[1]/$numamps;
  for( my $ct = 2; $ct < $repbins; ++$ct ) {
    printf "%6.2f%% amplicons have some barcodes at >~ %dx representation (%.2f%% at over 2*SD).$nlchr",
      100*$narcnts[$ct]/$numamps, (1<<($ct-1)), 100*$nar2sds[$ct]/$numamps;
  }
  print $nlchr;
}

# ------------- Generate the downsampling matrix -------------

# get total coverage and max duplication fold
my ($totcov,$maxdup) = (0,0);
for( my $an = 0; $an < $numamps; ++$an ) {
  my $row = $rows[$an];
  my $dsrow = $dsrows[$an];
  for( my $bc = 0; $bc < $numbcs; ++$bc ) {
    $totcov += $dsrow->[$bc];
    $maxdup = $row->[$bc] if( $row->[$bc] > $maxdup );
  }
}

# determine effective (1x) threshold for de-duping (if removing too many reads at 1x)
my $dupThres = 1;
unless( $dsonly || $dedupTarget <= 0 ) {
  my $targetreads = int(0.5+0.01*$dedupTarget*$totcov);
  $dupThres = getDedupThres($targetreads,$maxdup);
}

my $dstcov = 0;
my $fixcov = 0;
my $ovrcov = 0;
for( my $an = 0; $an < $numamps; ++$an ) {
  my $row = $rows[$an];
  my $dsrow = $dsrows[$an];
  my @ddcov = (0,0);
  for( my $bc = 0; $bc < $numbcs; ++$bc ) {
    my $origVal = $dsrow->[$bc];
    unless( $dsonly || $origVal <= $minDSCov ) {
      # do not down sample to median if already at or below that value
      $dsrow->[$bc] = int(0.5+$origVal*$dupThres/$row->[$bc]) if( $row->[$bc] > $dupThres );
      # allow sample to go down to but not below minimum sampling level
      $dsrow->[$bc] = $minDSCov if( $dsrow->[$bc] < $minDSCov );
    }
    # record number of reads at <= minimum coverage separate from those that might still be downsized
    $dstcov += $dsrow->[$bc];
    $fixcov += $dsrow->[$bc] if( $dsrow->[$bc] <= $minDSCov );
    $ovrcov += $dsrow->[$bc]-1 if( $dsrow->[$bc] > 1 );
    # restore original #reads for this barcode if pre-merged BC's option was employed
    if( $mergebcs ) {
      $ddcov[$bc & 1] += $dsrow->[$bc];
      $dsrow->[$bc] = $origVal;
    }
  }
  # perform random sampling for amplicon across merged barcodes (per strand)
  if( $mergebcs ) {
    for( my $dir = 0; $dir <= 1; ++$dir ) {
      my @deck;
      my $nbc = 0;
      for( my $bc = $dir; $bc < $numbcs; $bc += 2 ) {
        $nbc += $dsrow->[$bc];
        push( @deck, ($bc) x $dsrow->[$bc] );
      }
      my $nsam = $nbc - $ddcov[$dir]; # number to remove from all barcodes in this dir
      for( my $i = 0; $i < $nsam; ++$i ) {
        my $dsiz = $nbc-$i;
        my $samp = int(rand($dsiz));
        --$dsrow->[$deck[$samp]];
        $deck[$samp] = $deck[--$dsiz];
      }
    }
  }
}
if( $logstats ) {
  if( $dsonly ) {
    printf STDERR "Initial total reads assigned to targets: $totcov - min.cov. reads = $fixcov @ $minDSCov - 1.\n",
  } elsif( $logds ) {
    printf "Down-sample to %.4fx duplication threshold: $totcov -> $dstcov (%.2f%%)${nlchr}", $dupThres, 100*$dstcov/$totcov;
  } else {
    printf STDERR "Down-sample to ${dupThres}x duplication threshold: $totcov -> $dstcov (%.4f%%) - min.cov. reads = $fixcov @ $minDSCov - 1.\n",
      100*$dstcov/$totcov;
  }
  printf "Assumed pre-merged barcodes (initial coverage thresholds ignorred)${nlchr}" if( $mergebcs );
}

# reduce minimum down-sampled coverage to 1 (for any barcode for any strand of an amplicon)
$minDSCov = 1;

# try to linearly downsample over the variable (non-minimum) reads if necessary
# several tries are made since round-off means very hard to get exact number of reads required
$targetDSFrac *= 0.01; # convert % to fraction
my $dstarget = int(0.5+$targetDSFrac*$totcov);
my $doneds = $dstarget < $dstcov;
for( my $ntry = 1; $ntry <= $num_dstries && $dstarget < $dstcov; ++$ntry ) {
  # no minimum coverage when assuming pre-merged barcodes
  $ovrcov = $dstcov if( $mergebcs );
  # fraction of the re-scalable coverage to downsample: (target reads) = frc * (reads over-minimum) + (reads fixed at <= minimum)
  my $frc = ($dstarget + $ovrcov - $dstcov) / $ovrcov;
  my $last_dstcov = $dstcov;
  $dstcov = 0;
  $ovrcov = 0;
  for( my $an = 0; $an < $numamps; ++$an ) {
    my $row = $rows[$an];
    my $dsrow = $dsrows[$an];
    if( $mergebcs ) {
      # decrease sampling across all barcodes randomly to target level as if they were merged
      # make two passes to deal with fwd and rev reads separately
      for( my $dir = 0; $dir <= 1; ++$dir ) {
        my @deck;
        my $nbc = 0;
        for( my $bc = $dir; $bc < $numbcs; $bc += 2 ) {
          $nbc += $dsrow->[$bc];
          push( @deck, ($bc) x $dsrow->[$bc] );
        }
        my $nsam = $nbc - int(0.5+$frc*$nbc); # number to remove to get to desired level
        for( my $i = 0; $i < $nsam; ++$i ) {
          my $dsiz = $nbc-$i;
          my $samp = int(rand($dsiz));
          --$dsrow->[$deck[$samp]];
          $deck[$samp] = $deck[--$dsiz];
        }
      }
      for( my $bc = 0; $bc < $numbcs; ++$bc ) {
        $dstcov += $dsrow->[$bc];
      }
      next;
    }
    # since ratio only applied to barcodes at above $minDSCov
    # individual amplicons may not be downsampled to overall target but overall the total coverage will be
    for( my $bc = 0; $bc < $numbcs; ++$bc ) {
      if( $dsrow->[$bc] > $minDSCov ) {
        # coverage below original 1x representation may now be downsampled
        $dsrow->[$bc] = $minDSCov + int(0.5+$frc*($dsrow->[$bc]-$minDSCov));
      }
      $dstcov += $dsrow->[$bc];
      $ovrcov += $dsrow->[$bc]-$minDSCov if( $dsrow->[$bc] > $minDSCov );
    }
  }
  printf STDERR "Down-sample to target coverage #$ntry: $totcov -> $dstcov (%.4f%%).\n", 100*$dstcov/$totcov if( $logopt );
  last if( $dstcov == $last_dstcov );
}

# output downsampling matrix
for( my $an = 0; $an < $numamps; ++$an ) {
  my $row = $dsrows[$an];
  printf DSMATOUT "%s\t%s", $ampid[$an], $locus[$an];
  for( my $bc = 0; $bc < $numbcs; ++$bc ) {
    printf DSMATOUT "\t%g", $row->[$bc];
  }
  print DSMATOUT "\n";
}
close(DSMATOUT);

if( $logstats && $doneds ) {
  if( $logds ) {
    printf "Down-sampled amplicon reads %s per strand: $totcov -> $dstcov (%.2f%%).${nlchr}", 
      ($mergebcs ? "over all barcodes" : "to minimum 1 read per barcode"), 100*$dstcov/$totcov;
  } else {
    printf STDERR "Down-sample to target coverage: $totcov -> $dstcov (%.4f%%).\n", 100*$dstcov/$totcov;
  }
}

#-------------------------------------

sub median
{
  my @data = sort {$a <=> $b} @{$_[0]};
  my $ndat = scalar(@data);
  return if( $ndat <= 0 );
  return $data[0] if( $ndat == 1 );
  my $mdat = $ndat >> 1;
  return ($ndat & 1) ? $data[$mdat] : 0.5*($data[$mdat-1]+$data[$mdat]);
}

sub quartiles
{
  my @data = sort {$a <=> $b} @{$_[0]};
  my $ndat = scalar(@data);
  return () if( $ndat <= 0 );
  return ($data[0],$data[0],$data[0]) if( $ndat == 1 );
  my $odd = $ndat & 1;
  my $mdat = $ndat >> 1;
  my $md = $odd ? $data[$mdat] : 0.5*($data[$mdat-1]+$data[$mdat]);
  $ndat = ($quartile_method == 2) ? $mdat+$odd : $mdat;
  my $qdat = $ndat >> 1;
  my $q1 = $ndat & 1 ? $data[$qdat] : 0.5*($data[$qdat-1]+$data[$qdat]);
  $qdat += ($quartile_method == 2) ? $mdat : $mdat+$odd;
  my $q3 = $ndat & 1 ? $data[$qdat] : 0.5*($data[$qdat-1]+$data[$qdat]);
  return ($q1,$md,$q3);
}

sub avedev
{
  my $data = $_[0];
  my $ndat = scalar(@$data);
  return () if( $ndat <= 0 );
  return ($data->[0],0) if( $ndat == 1 );
  my $av = 0;
  for( my $i = 0; $i < $ndat; ++$i ) {
    $av += $data->[$i];
  }
  $av /= $ndat;
  my $var = 0;
  my $dif;
  for( my $i = 0; $i < $ndat; ++$i ) {
    $dif = $data->[$i] - $av;
    $var += $dif * $dif;
  }
  return ($av,sqrt($var/($ndat-1)) );
}

sub setRandSeed {
  my $sdstr = $_[0];
  my $seed = 0;
  for( my $i = 0; $i < length($sdstr); ++$i ) {
    $seed += ord( substr($sdstr,$i,1) );
  }
  return srand($seed);
}

sub getDedupThres {
  # binary search for target threshold approximately giving target reads
  my $trgrds = $_[0];
  my $maxdup = $_[1];
  my $maxacc = $_[2] ? $_[2] : 16;
  print STDERR "Max. dup = $maxdup\n" if( $logopt );
  # add extra levels of precision for high maxdup
  $maxacc += int($maxdup/64);
  my ($tlow,$thi) = (1,$maxdup);
  # first check 1x dedup level
  my $testrds = getToDedupThres(1);
  return 1 if( $testrds >= $trgrds );
  # binary search between low and high value until hit thereshold or too many cycles
  my ($thres,$ltestrds) = (1,0);
  while( --$maxacc >= 0 )
  {
    $thres = ($tlow + $thi)/2;
    $ltestrds = $testrds;
    $testrds = getToDedupThres($thres);
    print STDERR "getToDedupThres($thres) = $testrds vs. $trgrds\n" if( $logopt );
    print STDERR "getToDedupThres($thres) = $testrds vs. $trgrds\n";
    last if( $testrds == $trgrds || $ltestrds == $testrds );
    $testrds < $trgrds ? $tlow = $thres : $thi = $thres;
  }
  return $thres;
}

sub getToDedupThres {
  # simulates the total coverage count at a given fold threshold
  my $thr = $_[0];
  my $num = 0;
  for( my $an = 0; $an < $numamps; ++$an ) {
    my $row = $rows[$an];
    my $dsrow = $dsrows[$an];
    for( my $bc = 0; $bc < $numbcs; ++$bc ) {
      my $cnt = $dsrow->[$bc];
      if( $cnt > $minDSCov ) {
        $cnt = int(0.5+$cnt*$thr/$row->[$bc]) if( $row->[$bc] > $thr );
        $cnt = $minDSCov if( $cnt < $minDSCov );
      }
      $num += $cnt;
    }
  }
  return $num;
}

