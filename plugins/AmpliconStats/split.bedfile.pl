#!/usr/bin/perl -s

use strict;
use warnings;

use POSIX qw(strftime);

my $bedfile = $ARGV[0];
my $maxReadLength = $ARGV[1];
my $minRegionNumber = $ARGV[2];
my $maxSpreadLength = $ARGV[3];
my $splitName = $ARGV[4];
my $bamFile = $ARGV[5] || die "Usage: bedfile maxReadLength minRegionNumber maxSpreadLength splitName bamFile\n";

my @lines;
#read bed file to array of lines
open (BF, $bedfile) || die "Can't open $bedfile: $!\n";
while (<BF>) {
    next if (/^track/);
    push @lines, $_;
}
close (BF);

my $lct = 0;
my $veryStart = 0;
my $lastStart = 0;
my $lastEnd = 0;
my $lastChr = '';
my $currChr = 0;
my $currStart = 0;
my $currEnd = '';
my $splitCt = 0;
my %mergedRegions;

for my $i (0 .. $#lines) {
    my $line = $lines[$i];
    my @temp = split("\t", $line);
    $currChr = $temp[0];
    $currStart = $temp[1];
    $currEnd = $temp[2];

    if ($currChr eq $lastChr) { #same chromosome situations
	#two regions are <= maxReadLength, eg 1000, then merge into current mergedRegion
	if ($currStart - $lastEnd <= $maxReadLength) {
	} else {
	    #if it's maxSpreadLength apart, put into a new region
	    if ($currStart - $lastEnd > $maxSpreadLength) {
		$splitCt++;
	    } else {
		if ($#{$mergedRegions{$splitCt}} < $minRegionNumber) {
		} else {
		    $splitCt++; #why specify minRegionNumber?
		}
	    }
	}
    } else { #different chromosomes guarantee a new mergedRegion 
	$splitCt++;
    }
    push @{ $mergedRegions{$splitCt} }, $line;
    &recordLastCoordinates;
}

#clean merged regions, based on 1, each merged region should have minRegionNumber regions; 2, regions within a merged region should be <= maxSpreadLength apart

=hold
foreach my $k (sort {$a <=> $b} keys %mergedRegions) {
    my @temp = split("\t", $mergedRegions{$k}[0]);
    my $chrCoors = $temp[0] . ':' . $temp[1] . '-';
    @temp = split("\t", $mergedRegions{$k}[-1]);
    $chrCoors .= $temp[2];

    #my $date = strftime "%m/%d/%Y", localtime;
    #print $date;
    my $cmd = "samtools view -F 0x4 $bamFile $chrCoors > /dev/null";
    print "At ", `date`, "begins to run $cmd\n";
    #&runAndCheckCommand(0, $cmd);
}
=cut

#=hold
foreach my $k (keys %mergedRegions) {
    open (F, ">split.$k.bed") || die "Can't open split.$k.bed: $!\n";
    for my $m (0 .. $#{ $mergedRegions{$k} }) {
	print F $mergedRegions{$k}[$m];
    }
    close (F) || die "Can't close split.$k.bed: $!\n";
}
#=cut

sub recordLastCoordinates {
    $lastChr = $currChr;
    $lastStart = $currStart;
    $lastEnd = $currEnd;    
}

sub runAndCheckCommand {
    my ($con, $command) = @_;
    my $retVal = 0;
    #print($command, "\n");
    $retVal += system($command);
    print "At ", `date`, "return value now is $retVal\n";
    if ($retVal != 0) {
	if ($con) {
	    $retVal = 0;
	    print "Program encountered problem after running $command, but keep going.\n";
	} else {
	    print "Program stop after running $command\n";
	    exit(1);
	}
    }
}
