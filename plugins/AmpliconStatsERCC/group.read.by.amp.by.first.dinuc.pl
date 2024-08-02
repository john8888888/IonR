#!/usr/bin/perl -w

use strict;
use warnings;

use List::Util qw( min max ); #array min/max

my $csvf = $ARGV[0];
my $bcf = $ARGV[1] || die "Usage: example_ampliconSummary.csv barcodeList.txt\n";

#use JSON qw(decoded_json); #plan to parse results.json

my $dir = ".";
my @bc;
=hold
for my $i (1 .. 9) {
    my $b = 'IonXpress_00' . $i;
    push @bc, $b;
}
for my $i (10 .. 96) {
    my $b = 'IonXpress_0' . $i;
    push @bc, $b;
}
for my $i (1 .. 9) {
    my $b = 'IonXpress2_00' . $i;
    push @bc, $b;
}
for my $i (10 .. 99) {
    my $b = 'IonXpress2_0' . $i;
    push @bc, $b;
}
for my $i (0 .. 8) {
    my $b = 'IonXpress2_10' . $i;
    push @bc, $b;
}
=cut
open (BCF, $bcf) || die "Can't open $bcf: $!\n";
while (<BCF>) {
    next if !/^barcode/i;
    my @work = split(/,/, $_);
    push @bc, $work[1];
}
close (BCF);


my @alldinucs = qw(AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT);

my $primerFile = '/home/txu/projects/20130924.strand.bias.analysis/Exome_Bioinfo_Pool_1.tab';
my %insseq; #the fist dinuc is the wanted
my %dinucF;
my %ampFDinuc;
open (PF, $primerFile) || die "Can't open $primerFile: $!\n";
<PF>;
while (<PF>) {
    chomp;
    my @w = split(/\t/, $_);
    $insseq{$w[0]} = $w[22];
}
close (PF);

my %amp1;
my $totalAmpCt = 0;
#my $f = '/home/txu/projects/20130924.strand.bias.analysis/method.HiFi.60/Auto_user_P01-176-EqHiFi60_1_MDA_4183_9114.IonXpress_008.AmpStats_ampliconSummary.csv';
#should be the same, no difference
#my $f = '15337.IonXpress2_008.AmpStats_ampliconSummary.csv';
my $f = $csvf;
open (F, $f) || die "Can't open $f: $!\n";
<F>;
while (<F>) {
    chomp;
    my @w = split(/,/, $_);
    if ($insseq{$w[0]}) {
	my $fdinuc = substr($insseq{$w[0]}, 0, 2);
	$ampFDinuc{$fdinuc}{$w[0]} = 1;
	$dinucF{$w[0]} = $fdinuc;
	$totalAmpCt++;
	$amp1{$w[0]}{fdinuc} = $fdinuc;
    } else {
	print STDERR "no insert found for $w[0].\n";
    }
}
close (F);

my %expt;
foreach my $dc (sort {$a cmp $b} keys %ampFDinuc) {
    my @w = keys %{ $ampFDinuc{$dc} };
    $expt{$dc} = (scalar @w) * 100 / $totalAmpCt;
}
foreach my $dc (sort {$a cmp $b} keys %expt) {
    #print STDERR "$dc\t";
    #printf STDERR "%.4f\n", $expt{$dc};
}



my %table;
print "Dinuc\tExpted";
foreach my $thebc (sort {$a cmp $b} @bc) {
    print "\t$thebc";
    opendir (D, $dir);
    #group by amp
    #my @csv = grep { /MDA_4183_9114/ and /AmpStats_ampliconSummary.csv/ } readdir D;
    #group by barcode

    #my @csv = grep { /IonXpress2_$thebc/ and /AmpStats_ampliconSummary.csv/ } readdir D;
    my @csv = grep { /$thebc/ and /AmpStats_ampliconSummary.csv/ } readdir D;
    my %amp;
    print STDERR (scalar @csv), " files for $thebc combined.\n";
foreach my $f (@csv) {
    #print STDERR "working on $f\n";
    open (F, $f) || die "Can't open $f: $!\n";
    <F>;
    while (<F>) {
	chomp;
	my @w = split(/,/, $_);
	$amp{$w[0]}{frn} += $w[1];
	$amp{$w[0]}{rrn} += $w[2];
	$amp{$w[0]}{trn} += $w[3];
	$amp{$w[0]}{ffil} += $w[6];
	$amp{$w[0]}{fee} += $w[12];
    }
    close (F);
}
close (D);

my %dinucCount;
my $totalRead;
foreach my $dc (sort {$a cmp $b} @alldinucs) {
    foreach my $k (keys %amp) {
	if ($dc eq $amp1{$k}{fdinuc}) {
	    #$dinucCount{$dc} += $amp{$k}{frn};
	    #$totalRead += $amp{$k}{frn};
	    #$dinucCount{$dc} += $amp{$k}{ffil};
	    #$totalRead += $amp{$k}{ffil};
	    $dinucCount{$dc} += $amp{$k}{fee};
	    $totalRead += $amp{$k}{fee};
	}
    }
}

#print "group read counts by dinuc:\n";
foreach my $dc (keys %dinucCount) {
    #print "$dc\t";
    #printf "%.4f\t", $dinucCount{$dc} * 100 / $totalRead;
    #printf "%.4f", $expt{$dc};
    #print "\n";

    $table{$dc}{$thebc} = $dinucCount{$dc} * 100 / $totalRead;
}

}
print "\n";

foreach my $dc (keys %table) {
    print "$dc\t";
    printf "%.4f", $expt{$dc};
    foreach my $bc (sort {$a cmp $b} keys %{ $table{$dc} }) {
	printf "\t%.4f", $table{$dc}{$bc};
    }
    print "\n";
}

exit;

