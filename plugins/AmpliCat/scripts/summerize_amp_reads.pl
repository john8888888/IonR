#!/usr/bin/perl

#--------- Begin command arg parsing ---------

(my $CMD = $0) =~ s{^(.*/)+}{};
my $DESCR = "Create a list of statistics for all amplicons and all reads from the output of a
similar statistics for individual amplicons, as output by classify_amp_reads.pl run.";
my $USAGE = "Usage:\n\t$CMD [options] <amplicat file>";
my $OPTIONS = "Options:
  -h ? --help Display Help information.
  -L <N> Minimum amp_reads read value for considering an ampicon as not failing. Default: 20.
  -T <N> Total number of on-target reads. This is used to calculate the percercentage of assignable reads. Default: 0.
     If N <= 0 the sum of ovlp_reads is used, which may exceed the corect number if reads overlap multiple amplicons.";
  
my $lowreads = 20;
my $ontarget = 0;

my $help = (scalar(@ARGV) == 0);
while( scalar(@ARGV) > 0 )
{
  last if($ARGV[0] !~ /^-/);
  my $opt = shift;
  if($opt eq '-L') {$lowreads = shift;}
  elsif($opt eq '-T') {$ontarget = shift;}
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
elsif( $nargs != 1 )
{
  print STDERR "$CMD: Invalid number of arguments.";
  print STDERR "$USAGE\n";
  exit 1;
}

my $amplicatfile = shift;

#--------- End command arg parsing ---------

# Read assays: amplicons and primers

my ($sum_nreads,$sum_vreads,$sum_sreads,$sum_fpreads) = (0,0,0,0);
my ($failing,$outerprm,$highbias,$lowthru,$poorthru) = (0,0,0,0,0);
my ($nonu,$lowq,$fpps,$eepc,$pdpc) = (0,0,0,0,0);

open(AMPLICAT,$amplicatfile) || die "Cannot open input file $amplicatfile";
my $namps = -1;
while(<AMPLICAT>) {
  next if( ++$namps == 0 );
  chomp;
  my ($ampid,$locus,$novlps,$nreads,$vreads,$bias,$fcov,$rcov,$mq0,$mp10,$fpp,$eep,$pdp,$sreads) = split('\t',$_);
  $sum_nreads += $nreads;
  $sum_vreads += $vreads;
  $sum_sreads += $sreads;
  $sum_fpreads += int(0.5+0.01*$fpp*$vreads);
  ++$failing if( $vreads < $lowreads );
  ++$outerprm if( $nreads > 0 && $vreads/$nreads <= 0.9 );
  next if( $vreads <= 0 );
  $bias = 100 - $bias if( $bias > 50 );
  if( $vreads >= $lowreads ) {
    ++$highbias if( $bias < 20 );
    ++$lowthru if( $bias >= 30 && $fcov+$rcov < 120 );
    ++$poorthru if( $bias >= 40 && $fpp < 20 && $fcov+$rcov < 120 );
  }
  ++$nonu if( $mq0 >= 20 );
  ++$lowq if( $mp10 >= 20 );
  ++$eepc if( $eep >= 20 );
  ++$fpps if( $fpp >= 20 );
  ++$pdpc if( $pdp >= 20 );
}
close(AMPLICAT);

$sum_nreads = $ontarget if( $ontarget > 0 );

printf "Total assigned reads:       %10d\n", $sum_nreads;

# prevent div0 errors for no reads or amplicons
$sum_nreads = 1 if( $sum_nreads == 0 );
$sum_vreads = 1 if( $sum_vreads == 0 );
$namps = 1 if( $namps < 1 );

printf "Non-outer assigned reads:   %10d\n", $sum_vreads;

printf "Valid assigned reads:        %9.2f%%\n", 100*$sum_vreads/$sum_nreads;
printf "Short reads:                 %9.2f%%\n", 100*$sum_sreads/$sum_vreads;
printf "Premature attenuated reads:  %9.2f%%\n", 100*($sum_sreads-$sum_fpreads)/$sum_vreads;
printf "False primed reads:          %9.2f%%\n", 100*$sum_fpreads/$sum_vreads;
printf "Failing amplicons:           %9.2f%%\n", 100*$failing/$namps;
printf "Outer-primed amplicons:      %9.2f%%\n", 100*$outerprm/$namps;
printf "Non-uniquely read amplicons: %9.2f%%\n", 100*$nonu/$namps;
printf "Low-quality read amplicons:  %9.2f%%\n", 100*$lowq/$namps;
printf "High strand bias amplicons:  %9.2f%%\n", 100*$highbias/$namps;
printf "Low read-through amplicons:  %9.2f%%\n", 100*$lowthru/$namps;
printf "Poor read-through amplicons: %9.2f%%\n", 100*$poorthru/$namps;
printf "Short read amplicons:        %9.2f%%\n", 100*$eepc/$namps;
printf "False primed amplicons:      %9.2f%%\n", 100*$fpps/$namps;
printf "Primer-dimer amplicons:      %9.2f%%\n", 100*$pdpc/$namps;

