#!/usr/bin/perl

#--------- Begin command arg parsing ---------

(my $CMD = $0) =~ s{^(.*/)+}{};
my $DESCR = "Write a BED file from an existing well-formatted BED file that only contains non-overlapping targets.
See options below for more information. Output is to STDOUT.";
my $USAGE = "Usage:\n\t$CMD [options] <input BED file> <output BED file>";
my $OPTIONS = "Options:
  -h ? --help Display Help information
  -l Write progress Log information to STDERR.
  -n Write (only) the Number of targets in the output BED file to STDOUT.";

my $numberOut = 0;
my $logopt = 0;

my $help = (scalar(@ARGV) == 0);
while( scalar(@ARGV) > 0 )
{
  last if($ARGV[0] !~ /^-/);
  my $opt = shift;
  if($opt eq '-n') {$numberOut = 1;}
  elsif($opt eq '-l') {$logopt = 1;}
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

my $bedfile = shift;
my $outfile = shift;

#--------- End command arg parsing ---------

# First scaqn for overlapping targets

open(BEDFILE,$bedfile) || die "Cannot open targets file $bedfile";
my $num_targs = 0;
my ($chrid,$srt,$end,$ampid) = ("",0,0,"");
my (@ampchr,@ampsrt,@ampend,@ampid,@ampovlp);

while(<BEDFILE>) {
  next if( /^track / );
  ($chrid,$srt,$end,$ampid) = split('\t',$_);
  $ampchr[$num_targs] = $chrid;
  $ampsrt[$num_targs] = $srt+1;
  $ampend[$num_targs] = $end+0;
  $ampid[$num_targs] = $ampid;
  $ampovlp[$num_targs] = 0;
  ++$num_targs;
}
close(BEDFILE);
print STDERR "Loaded $num_targs targets\n" if( $logopt );

# identify overlapping targets
my $novlps = 0;
for( my $q = 0; $q < $num_targs-1; ++$q ) {
  next if( $ampovlp[$q] );
  my $qchr = $ampchr[$q];
  my $qsrt = $ampsrt[$q];
  my $qend = $ampend[$q];
  for( my $t = $q+1; $t < $num_targs; ++$t ) {
    next if( $ampovlp[$t] );
    next if( $ampchr[$q] ne $qchr );
    unless( $qsrt > $ampend[$t] || $qend < $ampsrt[$t] ) {
      printf STDERR "Overlap %d for $qchr:$qsrt-$qend and $ampchr[$q]:$ampsrt[$t]-$ampend[$t]\n",++$novlps if( $logopt );
      $ampovlp[$q] = 1;
      $ampovlp[$t] = 1;
      last;
    }
  }
}

# create list of non-overlapping targets

open(OUTFILE,">$outfile") || die "Cannot open targets output file $outfile";
open(BEDFILE,$bedfile) || die "Cannot open targets file $bedfile";
my ($num_in,$num_out) = (0,0);
while(<BEDFILE>) {
  if( /^track / ) {
    print OUTFILE;
    next;
  }
  unless( $ampovlp[$num_in] ) {
    print OUTFILE;
    ++$num_out;
  }
  ++$num_in;
}
close(OUTFILE);
close(BEDFILE);

print "$num_out\n" if( $numberOut );

