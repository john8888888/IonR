
#!/usr/bin/perl -w

#replacing shell script, in which 'cut -d"," -f1' does not work as expected

use strict;
use warnings;

my $ots = $ARGV[0]; #offtargetsummary
my $offbed = $ARGV[1]; #offtargetnamedbed
my $out = $ARGV[2] || die "Usage: input_off_target_summary_csv_file input_off_target_named_bed output_snap_tab_file\n";

my %off;
my %pos;
open (OB, $offbed) || die "Can't open $offbed: $!\n";
while (<OB>) {
    chomp;
    my @w = split(/\t/, $_);
    #$off{$w[3]} = "$w[4]\t$w[5]\t$w[6]\t$w[7]";
    $off{$w[3]} = "$w[4]\t$w[6]";
    $pos{$w[3]} = "$w[0]\t$w[1]\t$w[2]";
}
close (OB);

open (OUT, ">$out") || die "Can't open $out: $!\n";
#print OUT "AmpliconId\tNumFwd\tNumRev\tNumTotal\tNumSnapFwd\tNumSnapRev\tNumSnapTotal\tPercSnapFwd\tPercSnapRev\tPercSnapTotal\tCrossPrimerFwd\tFwdCIGAR\tCrossPrimerRev\tRevCIGAR\n";
#print OUT "AmpliconId\tChr\tStart\tEnd\tNumFwd\tNumRev\tNumTotal\tNumSnapFwd\tNumSnapRev\tNumSnapTotal\tPercSnapFwd\tPercSnapRev\tPercSnapTotal\tCrossPrimerFwd\tFwdCIGAR\tCrossPrimerRev\tRevCIGAR\n";
print OUT "AmpliconId\tChr\tStart\tEnd\tNumFwd\tNumRev\tNumTotal\tNumSnapFwd\tNumSnapRev\tNumSnapTotal\tPercSnapFwd\tPercSnapRev\tPercSnapTotal\tCrossPrimerFwd\tCrossPrimerRev\n";


open (OTS, $ots) || die "can't open $ots: $!\n";
while (<OTS>) {
    next if (/^AmpliconId/i);
    chomp;
    my @w = split(/,/, $_);
    my $ampid = $w[0];
    my $numfwd = $w[1];
    my $numrev = $w[2];
    my $numttl = $w[3];
    my $numsnapfwd = $w[22];
    my $numsnaprev = $w[23];
    my $numsnapttl = $numsnapfwd + $numsnaprev;
    my $cross = $off{$ampid};
    #print OUT "$ampid\t$numfwd\t$numrev\t$numttl\t$numsnapfwd\t$numsnaprev\t$numsnapttl";
    print OUT "$ampid\t$pos{$ampid}\t$numfwd\t$numrev\t$numttl\t$numsnapfwd\t$numsnaprev\t$numsnapttl";
    if ($numfwd == 0 || $numrev == 0) {
	print OUT "\tNA\tNA\tNA";
    } else {
	my $precsnapfwd = $numsnapfwd * 100 / $numfwd;
	my $precsnaprev = $numsnaprev * 100 / $numrev;
	my $precsnapttl = $numsnapttl * 100 / $numttl;
	printf OUT "\t%4.2f\t%4.2f\t%4.2f", $precsnapfwd, $precsnaprev, $precsnapttl;
    }
    print OUT "\t$cross\n";
}
close (OTS);
close (OUT);
