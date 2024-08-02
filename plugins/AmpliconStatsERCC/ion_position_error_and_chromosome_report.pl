#!/usr/bin/env perl

use strict;
use warnings;
use FileHandle;

my $bed = $ARGV[0];
my $bam = $ARGV[1];
#my $samtools = $ARGV[2] || die "Usage: $0 bed.file bam.file samtools\n";
my $coverage = $ARGV[2];
my $errorOut = $ARGV[3];
my $chrRepOut = $ARGV[4] || die "Usage: $0 bed.file bam.file coverage.file position.error.outfile chromosome.report.outfile\n";

my $total_base = 0;
my $total_err = 0;
my $total_del = 0;
my $total_ins = 0;
my $total_mis = 0;
my $accuracy = 0;
my $ct = 0;

my %chromosome;

my $sampipe = "samtools view -F 0x4 -L $bed $bam 2> /dev/null |";
open (IN, $sampipe) || die "Can't read samtools pipe: $!\n\t-- $sampipe\n";

while (<IN>) {
    next if (/^@/);
    $ct++;
    print STDERR "$ct\r" if ($ct % 10000 == 0);
    my @temp = split(/\t/, $_);
    next if ($temp[0] eq 'name');

    #chromosome
    if ($temp[1] == 0) {
	$chromosome{$temp[2]}{fwd}++;
    } elsif ($temp[1] == 16) {
	$chromosome{$temp[2]}{rev}++;
    } else {}


    #cigar
    my $match = $temp[5];
    my $match1 = $match;
    #$match1 =~ s/[MDIS]/:/g;
    $match1 =~ s/[MDISH]/:/g;
    my @dig = split(/:/, $match1);
    my @cha = split(/\d+/, $match);
    die "Not expected cigar string $match\n" if ($cha[0] ne '');
    for my $i (1 .. $#cha) {
	my $j = $i - 1;
	if ($cha[$i] eq 'M') {
	    $total_base += $dig[$j];
	} elsif ($cha[$i] eq 'I') {
	    $total_base += $dig[$j];
	    $total_ins += $dig[$j];
	} elsif ($cha[$i] eq 'D') {
	    $total_del += $dig[$j];
	} elsif ($cha[$i] eq 'S') {
	    #soft/hard clipping
	    #die "Not expected cigar character in $match because of $cha[$i] not MIDS\n";
	}
    }

    #MD field
    my $mdf = 0;
    for my $i (10..$#temp) {
        if( $temp[$i] =~ m/^MD:/ ) {
            $mdf = $i;
            last;
        }
    }
    if( $mdf > 0 )
    {
        my @work = split(/\:/, $temp[$mdf]);
        my $mdstr = $work[2]; #MD string: reg exp, [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
        my $mis = $mdstr;
        $mis =~ s/\^[ACGT]+//g;
        $mis =~ s/\d//g;
        $total_mis += length($mis);
        #if mismatches at the end
        my $endmis = $mdstr;
        if ($endmis =~ /([ACGT]0)+$/) {
            #print "$mdstr has end mis of $&\n";
            my $endmisstr = $&;
            $endmisstr =~ s/0//g;
            $total_mis -= length($endmisstr);
            $total_ins += length($endmisstr);
        }
    }
}
#close (IN)  || die "Can't close $in: $!\n";

$total_err = $total_ins + $total_del + $total_mis;
$accuracy = (1-($total_err/$total_base))*100;

open (EO, ">$errorOut") || die "Can't open $errorOut: $!\n";
print EO "Total bases called = $total_base\n";
print EO "Number of insertions = $total_ins\n";
print EO "Number of mismatches = $total_mis\n";
print EO  "Number of deletions = $total_del\n";
print EO "Total Errors = $total_err\n";
print EO "Average Read Accuracy = $accuracy\n";
close (EO);

#coverage data covers all bases, including base with zero coverage
open (C, $coverage) || die "Can't open $coverage: $!\n";
while (<C>) {
    chomp;
    my @w = split(/\t/, $_);
    push @{ $chromosome{$w[0]}{depth}}, $w[2];
}
close (C);

foreach my $c (keys %chromosome) {
    foreach my $d (@{ $chromosome{$c}{depth} }) {
	$chromosome{$c}{totalChrBase} += $d;
	$chromosome{$c}{'1x'}++ if ($d >= 1);
	$chromosome{$c}{'10x'}++ if ($d >= 10);
	$chromosome{$c}{'100x'}++ if ($d >= 100);
	$chromosome{$c}{'500x'}++ if ($d >= 500);
	$chromosome{$c}{'1000x'}++ if ($d >= 1000);
    }
}
open (CRO, ">$chrRepOut") || die "Can't open $chrRepOut: $!\n";
print CRO "Chrom,BasesCovered,TotalBases,FwdReads,RevReads,TotalReads,1x,10x,100x,500x,1000x\n";
foreach my $c (sort sortChr (keys %chromosome)) {
    #my $totalReads = $chromosome{$c}{fwd} + $chromosome{$c}{rev};
    my $totalReads = $chromosome{$c}{fwd} ? $chromosome{$c}{fwd} : 0 + $chromosome{$c}{rev} ? $chromosome{$c}{rev} : 0;
    my $baseCovered = scalar @{ $chromosome{$c}{depth} };
    my $c1x = $chromosome{$c}{'1x'} ? $chromosome{$c}{'1x'} : 0;
    my $c10x = $chromosome{$c}{'10x'} ? $chromosome{$c}{'10x'} : 0;
    my $c100x = $chromosome{$c}{'100x'} ? $chromosome{$c}{'100x'} : 0;
    my $c500x = $chromosome{$c}{'500x'} ? $chromosome{$c}{'500x'} : 0;
    my $c1000x = $chromosome{$c}{'1000x'} ? $chromosome{$c}{'1000x'} : 0;
    #printf CRO "%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", ($c, $baseCovered, $chromosome{$c}{totalChrBase}, $chromosome{$c}{fwd}, $chromosome{$c}{rev}, $totalReads, $c1x, $c10x, $c100x, $c500x, $c1000x);
    printf CRO "%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", ($c, $baseCovered, $chromosome{$c}{totalChrBase} ? $chromosome{$c}{totalChrBase} : 0, $chromosome{$c}{fwd} ? $chromosome{$c}{fwd} : 0, $chromosome{$c}{rev} ? $chromosome{$c}{rev} : 0, $totalReads, $c1x, $c10x, $c100x, $c500x, $c1000x);
}
close (CRO);


sub sortChr {
    my $number_a = substr($a, 3);
    $number_a =~ s/[A-Z]//g;
    my $number_b = substr($b, 3);
    $number_b =~ s/[A-Z]//g;
    my $letter_a = substr($a, 3);
    $letter_a =~ s/\d//g;
    my $letter_b = substr($b, 3);
    $letter_b =~ s/\d//g;

    #my $number_b = substr($b, 3) =~ s/[A-Z]//g;
    #my $letter_b = substr($b, 3) =~ s/\d//g;
    #print "$number_a, $number_b, $letter_a, $letter_b\n";
    return $number_a <=> $number_b or $letter_a cmp $letter_b;
}
