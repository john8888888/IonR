#!/usr/bin/env perl

use strict;
use warnings;
use FileHandle;

my $bed = $ARGV[0];
my $bam = $ARGV[1];
my $samtools = $ARGV[2] || die "Usage: $0 bed.file bam.file samtools\n";

my $total_base = 0;
my $total_err = 0;
my $total_del = 0;
my $total_ins = 0;
my $total_mis = 0;
my $accuracy = 0;
my $ct = 0;

my $sampipe = "$samtools view -F 0x4 -L $bed $bam 2> /dev/null |";
open (IN, $sampipe) || die "Can't read samtools pipe: $!\n\t-- $sampipe\n";

while (<IN>) {
    next if (/^@/);
    $ct++;
    print STDERR "$ct\r" if ($ct % 10000 == 0);
    my @temp = split(/\t/, $_);
    next if ($temp[0] eq 'name');

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
print "Total bases called = $total_base\n";
print "Number of insertions = $total_ins\n";
print "Number of mismatches = $total_mis\n";
print "Number of deletions = $total_del\n";
print "Total Errors = $total_err\n";
print "Average Read Accuracy = $accuracy\n";

