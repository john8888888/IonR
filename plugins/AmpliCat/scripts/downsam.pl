#!/usr/bin/perl
# Copyright (C) 2012 Ion Torrent Systems, Inc. All Rights Reserved

use File::Basename;

(my $CMD = $0) =~ s{^(.*/)+}{};
my $DESCR = "Create SAM file from a merged set of sampled barcode files on a per-strand per-amplicon basis.
By default, up to the best 2 reads (1 per strand) will be chosen as longest aligned start, with best quality and nearest
start to insert start to split ties. Or a barcode down-sampling file will be used if supplied.
The files provided are read-to-amplicon mapping files for each barcode, which is also the root the BAM files and must
start with the barcode ID for sake of matching to the barcode sampling defined by the down-sampling file, if provided.";

my $USAGE = "Usage:\n\t$CMD [options] [-D <down-sampling file>] <Out bam file name> <assigned read mapping file 1> [<file2> ...]\n";
my $OPTIONS = "Options:
  -h ? --help Display Help information
  -l Output extra run Log information to STDERR.
  -v Output Verbose run log information to STDERR.
  -b Use Best single alignment per amplicon per strand: Longest aligned then by average QV. (-D option becomes redundant)
  -B <file> BED file containing target locations. Used to prefer reads better aligned to inserts for AmpliSeq with -b option.
     BED file must include amplicon ID as 4th field.
  -D <file> Down-sampling matrix. Default: '' => Merge all barcode reads.
  -S <int>  Stop after writing results for this many amplicons. Use 0 to stop after writing header. Default: -1 => no stop.";

my $logopt = 0;
my $maxSample = -1;

my $bestreads = 0;
my $dsmatrixfile = '';
my $targetsfile = '';

my $maxAmpOvlp = 10;	# defines number of amplicons (reads) buffered before first on stack is output

my $help = (scalar(@ARGV) == 0);
while( scalar(@ARGV) > 0 )
{
  last if($ARGV[0] !~ /^-/);
  my $opt = shift;
  if($opt eq '-l') {$logopt = 1;}
  elsif($opt eq '-v') {$logopt = 2;}
  elsif($opt eq '-b') {$bestreads = 1;}
  elsif($opt eq '-B') {$targetsfile = shift;}
  elsif($opt eq '-D') {$dsmatrixfile = shift;}
  elsif($opt eq '-S') {$maxSample = int(shift);}
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
elsif( $nargs < 2 )
{
  print STDERR "Error: Insuficient args: Usage:\n  $CMD [options] [-D <down-sampling matrix file>] <Out bam file name> <assigned read mapping file 1> [<file2> ...]\n";
  exit 1;
}

my $outfileroot = shift;
$outfileroot =~ s/\.bam$//;
my $samout = $outfileroot.".sam";

$dsmatrixfile = '' if( $dsmatrixfile eq '-' || $bestreads );
my $haveDsmatrix = $dsmatrixfile ne '' ? 1 : 0;

$targetsfile = '' if( $targetsfile eq '-' );
my $haveTargets = $targetsfile ne '' ? 1 : 0;

# temporary performance enhancement (until this option is also employed for biased down-samping)
$haveTargets = 0 if( $haveDsmatrix );

my $cacheAllReads = ($haveDsmatrix && $bestreads == 0) ? 1 : 0;

print STDERR "$CMD running..\n" if( $logout );
printf STDERR "dsmatrixfile = '%s'\nsamout = '$samout'\nNum input files: %d\n",($bestreads ? "ignorred" : $dsmatrixfile),($nargs-1) if( $logopt > 1 );
print STDERR "bestreads = $bestreads, haveDsmatrix = $haveDsmatrix, haveTargets = $haveTargets, cacheAllReads = $cacheAllReads\n" if( $logopt > 1 );

# -------------- End command arg processing --------------

# read target insert ranges
my (%ampsrts,%ampends);
if( $haveTargets ) {
  open(BEDFILE,$targetsfile) || die "Could not open targets BED file '$targetsfile'.\n";
  my $namps = 0;
  while(<BEDFILE>) {
    chomp;
    my ($bchr,$bsrt,$bend,$bname) = split('\t',$_);
    next if( $bchr !~ /\S/ );
    next if( $bchr =~ /^track / );
    ++$namps;
    $ampsrts{$bname} = 1+$bsrt;
    $ampends{$bname} = 0+$bend;
  }
  close(BEDFILE);
  print STDERR "$CMD: Read $namps amplicons from targets BED file.\n" if( $logopt );
}

# Read entire down-sampling matrix into memory (needed for look up)
my (%ampfield,%ampdsns);
if( $haveDsmatrix ) {
  my $namps = -1;
  open(DSMATRIX,$dsmatrixfile) || die "Could not open down-sampling matrix file '$dsmatrixfile'.\n";
  while(<DSMATRIX>) {
    chomp;
    my @fields = split('\t',$_);
    if( ++$namps == 0 ) {
      for( my $i = 2; $i < scalar(@fields); ++$i ) {
        $ampfield{$fields[$i]} = $i-2;
      }
      next;
    }
    my $ampid = shift(@fields);
    my $locus = shift(@fields);
    $ampdsns{$ampid} = \@fields;
    # sanity check
    if( $haveTargets && !defined($ampsrts{$ampid}) ) {
      print STDERR "$CMD: Warning: Reads for amplicon $ampid but no insert defined in BED file supplied.\n";
    }
  }
  close(DSMATRIX);
  print STDERR "$CMD: Read $namps amplicons from down-sampling matrix file.\n" if( $logopt );
}

# make one pass through the list of (BAM) file to create a merged header file
my $rxsub = 's/(^\@RG\s.*)(\sSM:)(.*?)($|\t)/$1${2}None$4/';
for( my $argn = 0; $argn <= $#ARGV; ++$argn ) {
  my $bamfile = $ARGV[$argn];
  $bamfile =~ s/\.downsample.map$//;
  $bamfile .= ".bam";
  # copy header from just first file
  unless( $argn ) {
    if( system("samtools view -H \"$bamfile\" | perl -pe '$rxsub' > \"$samout\"") ) {
      print STDERR "Warning: Could not open bamfile '$bamfile' to create primary SAM header.\n";
      exit 1;
    }
    next;
  }
  # extract and modify just the RG lines from the other headers
  if( system("samtools view -H \"$bamfile\" | grep -e '^\@RG' | perl -pe '$rxsub' >> \"$samout\"") ) {
    print STDERR "Warning: Could not open bamfile '$bamfile' to add RG lines to header.\n";
  }
}
unless( $maxSample ) {
  print STDERR "End after writing merged headers for option -S 0\n";
  exit(0);
}
open(SAMOUT,">>$samout") || die "Could not open to append to merged SAM header file '$samout'.\n";

# amplicon read stack buffer - allows reads to be interlaced vs. amplicon assignments
my (@amprds_fwd,@amprds_rev,%amptable);
my @ampnames = ('') x $maxAmpOvlp;
my $ampslot = 0;
my (@bamp_len_fwd,@bamp_len_rev,@bamp_aqv_fwd,@bamp_aqv_rev);

# these varibles track current barcode and are also global to buffering routines
my ($barcode,$bcfwd,$bcrev,$bcfwdidx,$bcrevidx);

# high level stats
my ($filesread,$samRecsWritten,$totalSamRecsWritten) = (0,0,0);

while( scalar(@ARGV) > 0 ) {
  # get input mapping and bam files
  my $sammap = shift;
  unless( -e $sammap ) {
    print STDERR "Could not locate read-to-target mapping file '$sammap' - skipping.\n";
    next;
  }
  open(SAMMAP,$sammap) || die "Could not open SAM read-to-target mapping file '$sammap'\n";
  my $bamfile = $sammap;
  $bamfile =~ s/\.downsample.map$//;
  $bamfile .= ".bam";
  open(BAMREAD,"samtools view -F 4 \"$bamfile\" |") || die "Could not locate BAM file '$bamfile'\n";
  ++$filesread;
  # use input file names for random-seed setting - hence results should be the same if same ban files re-analyzed
  $bamfile =~ s/.*\///;
  setRandSeed($CMD.$bamfile);
  # get barcode to find indexes into amplicon sample numbers
  $bamfile =~ m/^(.*?_\d+?)_/;
  $barcode = $1;
  print STDERR "Processing reads for barcode $barcode...\n" if( $logopt > 1 );
  if( $haveDsmatrix ) {
    $bcfwd = $barcode."_fwd";
    $bcrev = $barcode."_rev";
    unless( defined($ampfield{$bcfwd}) ) {
      print STDERR "$CMD: Could not locate read-to-target mapping data for barcode key $bcfwd.\n";
      print STDERR "  - skipping file $bamfile\n";
      close(BAMREAD);
      close(SAMMAP);
      next;
    }
    $bcfwdidx = $ampfield{$bcfwd};
    $bcrevidx = $ampfield{$bcrev};
  }
  # extract reads listed in the read-to-amplicon mapping file
  my $line;
  my $nrids = 0;
  while(<SAMMAP>) {
    chomp;
    my ($rid,$aid) = split('\t');
    # BAM reads should be in same ordering as mapping file
    my $dir = 0;
    while(<BAMREAD>) {
      chomp;
      $line = $_;
      /^(.+?)\t(\d+?)\t/;
      if( $1 eq $rid ) {
        $dir = ($2 & 16) ? -1 : 1;
        last;
      }
    }
    die "Read-to-target mapping file out of alignment with BAM file '$bamfile' for read '$rid' after reading line:\n$lastLine\n" unless( $dir );
    # collect all reads per amplicon and output when sure not more reads out of phase with amplicons
    bankRead($aid,$dir,$line);
  }
  # write the reads for last amplicons
  flushBank();
  close(BAMREAD);
  close(SAMPMAP);
  print STDERR "Wrote $samRecsWritten reads for barcode $barcode.\n" if( $logopt );
  $totalSamRecsWritten += $samRecsWritten;
  $samRecsWritten = 0;
}
close(SAMOUT);
print STDERR "Wrote total of $totalSamRecsWritten for $filesread input files\n" if( $logopt );

# ------------------------- END -----------------------

sub setRandSeed {
  my $sdstr = $_[0];
  my $seed = 0;
  for( my $i = 0; $i < length($sdstr); ++$i ) {
    $seed += ord( substr($sdstr,$i,1) );
  }
  return srand($seed);
}

sub writeSampleAmpReads {
  my $list = $_[0];
  my $nsam = $_[1];
  my $nrds = scalar(@$list);
  # do not randomize samping if using all (merged file still needs sorting)
  if( $nsam >= $nrds ) {
    for( my $i = 0; $i < $nrds; ++$i ) {
      print SAMOUT "$list->[$i]\n";
    }
    $samRecsWritten += $nrds;
    printf STDERR "Wrote all %d reads for $_[2]. Total: $samRecsWritten.\n",$nrds if( $logopt > 1 );
  } else {
    my @deck = (0..$nrds);
    $nsam = $nrds if( $nrds < $nsam );
    for( my $i = 0; $i < $nsam; ++$i ) {
      my $dsiz = $nrds-$i;
      my $samp = int(rand($dsiz));
      print SAMOUT "$list->[$deck[$samp]]\n";
      $deck[$samp] = $deck[--$dsiz];
    }
    $samRecsWritten += $nsam;
    print STDERR "Wrote $nsam reads of $nrds for $_[2]. Total: $samRecsWritten.\n" if( $logopt > 1 );
  }
  if( $maxSample > 0 && $samRecsWritten >= $maxSample ) {
    print STDERR "Stopped after $samRecsWritten reads written for $filesread input files due to option -S $maxSample.\n";
    exit(0);
  }
}

sub bankRead {
  my $ampid = $_[0];
  my $rddir = $_[1];
  my $samrc = $_[2];
  # check if this amplicon is curently in buffered list
  my $bank;
  unless( defined($amptable{$ampid}) ) {
    # rotate next bank pointer
    $ampslot = 0 if( ++$ampslot >= $maxAmpOvlp );
    # if bank box used up then sample and write
    my $oaid = $ampnames[$ampslot];
    if( $oaid ne '' ) {
      if( $haveDsmatrix ) {
        unless( defined($ampdsns{$oaid}) ) {
          print STDERR "$CMD: Could not locate down-sample data for amplicon $ampid\n";
          return;
        }
        writeSampleAmpReads($amprds_fwd[$ampslot],$ampdsns{$oaid}->[$bcfwdidx],$oaid);
        writeSampleAmpReads($amprds_rev[$ampslot],$ampdsns{$oaid}->[$bcrevidx],$oaid);
      } else {
        writeSampleAmpReads($amprds_fwd[$ampslot],1,$oaid);
        writeSampleAmpReads($amprds_rev[$ampslot],1,$oaid);
        $bamp_len_fwd[$ampslot] = $bamp_len_rev[$ampslot] = 0;
        $bamp_aqv_fwd[$ampslot] = $bamp_aqv_rev[$ampslot] = 0;
      }
      $amprds_fwd[$ampslot] = [];
      $amprds_rev[$ampslot] = [];
      delete($amptable{$oaid});
    }
    $amptable{$ampid} = $ampslot;
    $ampnames[$ampslot] = $ampid;
    $bank = $ampslot;
  } else {
    $bank = $amptable{$ampid};
  }
  if( $cacheAllReads ) {
    $rddir > 0 ? push( @{$amprds_fwd[$bank]}, $samrc ) : push( @{$amprds_rev[$bank]}, $samrc );
  } else {
    # just keep best read per strand - here get aligned length (for position and selector)
    my @fields = split("\t",$samrc);
    my $cig = $fields[5];
    my $len = 0;
    while( $cig =~ s/^(\d+)(.)// ) {
      $len += $1 if( $2 eq "M" || $2 eq "D" || $2 eq "X" || $2 eq "=" );
    }
    # modify effective read length if insert coords are available
    if( $haveTargets && defined($ampsrts{$ampid}) ) {
      my $cor = $fields[3] - $ampsrts{$ampid};
      $len += $cor if( $cor < 0 );
      $cor = $ampends{$ampid} - ($fields[3]+$len+1);
      $len += $cor if( $cor < 0 );
    }
    # primary selector is read length, after clipping length to inserts
    my $bestLen = $rddir > 0 ? $bamp_len_fwd[$bank] : $bamp_len_rev[$bank];
    unless( $len < $bestLen ) {
      # split ties by average QV => also need to save for best length
      # - average QV includes soft-clips and INDELs but should not matter for longest alignments
      my @qva = unpack("C*",$fields[10]);
      my $aqv = 0;
      for( my $i = 0; $i < scalar(@qva); ++$i ) { $aqv += $qva[$i] }
      $aqv /= scalar(@qva);
      my $bestAqv = $rddir > 0 ? $bamp_aqv_fwd[$bank] : $bamp_aqv_rev[$bank];
      if( $len > $bestlen || $qvs > $bestAqv ) {
        if( $rddir > 0 ) {
          $amprds_fwd[$bank][0] = $samrc;
          $bamp_len_fwd[$bank] = $len;
          $bamp_aqv_fwd[$bank] = $aqv;
        } else {
          $amprds_rev[$bank][0] = $samrc;
          $bamp_len_rev[$bank] = $len;
          $bamp_aqv_rev[$bank] = $aqv;
        }
      }
    }
  }
}

sub flushBank {
  # dump out in order stored
  my $bank = $ampslot;
  do {
    $ampslot = 0 if( ++$ampslot >= $maxAmpOvlp );
    my $oaid = $ampnames[$ampslot];
    if( $oaid ne '' ) {
      if( $haveDsmatrix ) {
        writeSampleAmpReads($amprds_fwd[$ampslot],$ampdsns{$oaid}->[$bcfwdidx],$oaid);
        writeSampleAmpReads($amprds_rev[$ampslot],$ampdsns{$oaid}->[$bcrevidx],$oaid);
      } else {
        writeSampleAmpReads($amprds_fwd[$ampslot],1,$oaid);
        writeSampleAmpReads($amprds_rev[$ampslot],1,$oaid);
        $bamp_len_fwd[$ampslot] = $bamp_len_rev[$ampslot] = 0;
        $bamp_aqv_fwd[$ampslot] = $bamp_aqv_rev[$ampslot] = 0;
      }
      $amprds_fwd[$ampslot] = [];
      $amprds_rev[$ampslot] = [];
    }
  } while( $ampslot != $bank );
  %amptable = ();
  @ampnames = ();
  $ampslot = 0;
}

