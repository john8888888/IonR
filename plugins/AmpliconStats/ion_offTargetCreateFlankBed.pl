#!/usr/bin/perl -w

#to replace the part of reading input bedfile, which could be large, in script offTarget.sh

use strict;
use warnings;
use POSIX qw(strftime);

my $bamFile = $ARGV[0];
my $bedFile = $ARGV[1];
my $outputDir = $ARGV[2];
my $baseDir = $ARGV[3];
my $reference = $ARGV[4];
my $analysis = $ARGV[5];
my $ampSumCSV = $ARGV[6];

my $date = localtime;
print "Beginning ion_offTargetCreateFlankBed.pl at $date\n";

my $gapAllowed = 10;
my $primerLen = 25;
my $primerInt = 10;

my $offTargetName = "$outputDir/$analysis.offtarget";
my $offTargetBam = "$offTargetName.bam";
my $offTargetBed = "$offTargetName.bed";
my $offTargetFa = "$offTargetName.fa";
my $offTargetBedNamed = "$offTargetName.named.bed";
my $proximityBed = "$offTargetName.proximity.bed";
my $offTargetAmpReport = "$offTargetName.ampSummary.csv";
my $offTargetSnap = "$offTargetName.snap.tab";
my $targetName = "$outputDir/$analysis.target";
my $targetFlankBed = "$targetName.flanks.bed";
my $offTargetFlankBed = "$offTargetName.flanks.bed";
my $targetFlankFa = "$targetName.flanks.fa";
my $offTargetFlankFa = "$offTargetName.flanks.fa";
my $offTargetAlign = "$offTargetName.flanks.align.sam";
my $offTargetFlankBedAligned = "$offTargetName.flanks.align.bed";



open (TARGETFLANKBED, ">$targetFlankBed") || die "Can't open $targetFlankBed: $!\n";
open (OFFTARGETFLANKBED, ">$offTargetFlankBed") || die "Can't open $offTargetFlankBed: $!\n";
open (OFFTARGETFLANKBEDALIGNED, ">$offTargetFlankBedAligned") || die "Can't open $offTargetFlankBedAligned: $!\n";
open (OFFTARGETBEDNAMED, ">$offTargetBedNamed") || die "Can't open $offTargetBedNamed: $!\n";


my $primerSeqCheck = 0;
my %ampSum;
my @definitions;
my $fwdPrimeFieldNum = 0;
my $revPrimeFieldNum = 0;
open (AS, $ampSumCSV) || die "Can't open $ampSumCSV: $!\n";
while (<AS>) {
    chomp;
    if (/^AmpliconId,NumFwd,NumRev/) {
	@definitions = split(/,/, $_);
	for my $i (0 .. $#definitions) {
	    $primerSeqCheck = 1 if (uc($definitions[$i]) eq "FWD_UPRIMER");
	    $fwdPrimeFieldNum = $i if (uc($definitions[$i]) eq "FWD_UPRIMER");
	    $revPrimeFieldNum = $i if (uc($definitions[$i]) eq "REV_UPRIMER");
	}
    } else {
	my @w = split(/,/, $_);
	if ($primerSeqCheck > 0) {
	    $ampSum{$w[0]}{fwdplen} = length($w[$fwdPrimeFieldNum]);
	    $ampSum{$w[0]}{revplen} = length($w[$revPrimeFieldNum]);
	} else {
	    $ampSum{$w[0]}{fwdplen} = $primerLen;
	    $ampSum{$w[0]}{revplen} = $primerLen;
	}
    }
}
close (AS);

$date = localtime;
print "Finished reading in ampSumCSV in ion_offTargetCreateFlankBed.pl at $date\n";

open (BF, $bedFile) || die "Can't open $bedFile: $!\n";
while (<BF>) {
    next if (/track/);
    chomp;
    my @w = split(/\t/, $_);
    my $chr = $w[0];
    my $start = $w[1];
    my $stop = $w[2];
    my $ampid = $w[3];
    my $newStart = $start - $ampSum{$ampid}{fwdplen}; #fwdplen may be primerLen
    my $newStop = $stop + $ampSum{$ampid}{revplen}; #revplen may be primerLen
    my $intStart = $start + $primerInt;
    my $intStop = $stop - $primerInt;
    
    print TARGETFLANKBED "$chr\t$newStart\t$intStart\t$ampid", "_fwd\n";
    print TARGETFLANKBED "$chr\t$intStop\t$newStop\t$ampid", "_rev\n";
}
close (BF);
close (TARGETFLANKBED);

$date = localtime;
print "Finished reading in bedFile in ion_offTargetCreateFlankBed.pl at $date\n";
print "Finished creating targetFlankBed in ion_offTargetCreateFlankBed.pl at $date\n";

my $count = 0;
open (OFFTARGETBED, $offTargetBed) || die "Can't open $offTargetBed: $!\n";
while (<OFFTARGETBED>) {
    $count++;
    chomp;
    my @w = split(/\t/, $_);
    my $chr = $w[0];
    my $start = $w[1];
    my $stop = $w[2];
    my $newStart = $start - $primerLen;
    my $intStart = $start + $primerLen;
    my $newStop = $stop + $primerLen;
    my $intStop = $stop - $primerLen;
    print OFFTARGETFLANKBED "$chr\t$newStart\t$intStart\t$count", "_fwd\n";
    print OFFTARGETFLANKBED "$chr\t$intStop\t$newStop\t$count", "_rev\n";
}
close (OFFTARGETBED);
close (OFFTARGETFLANKBED);

$date = localtime;
print "Finished creating offTargetFlankBed in ion_offTargetCreateFlankBed.pl at $date\n";

my $cmd = "$baseDir/fastaFromBed -name -fi $reference -bed $targetFlankBed -fo $targetFlankFa";
runAndCheckCommand(1);
$cmd = "$baseDir/fastaFromBed -name -fi $reference -bed $offTargetFlankBed -fo $offTargetFlankFa";
runAndCheckCommand(1);

$cmd = "tmap index -f $targetFlankFa";
runAndCheckCommand(1);
$cmd = "tmap mapall -f $targetFlankFa -r $offTargetFlankFa -a 2 -g 0 stage1 map1 map2 map3 | $baseDir/samtools view -S - -F 0x4 > $offTargetAlign";
runAndCheckCommand(1);

$date = localtime;
print "Finished tmap-ping off target reads in ion_offTargetCreateFlankBed.pl at $date\n";

my %alignment;
open (OTA, $offTargetAlign) || die "Can't open $offTargetAlign: $!\n";
while (<OTA>) {
    chomp;
    my @w = split(/\t/, $_);
    push @{ $alignment{$w[0]}{hits} }, $w[2]; #different than offTarget.sh
    $alignment{$w[0]}{cigar} = $w[5];
}
close (OTA);

$date = localtime;
print "Finished reading in tmap alignments in ion_offTargetCreateFlankBed.pl at $date\n";

$count = 0;
my %entries; #possibly multiple entries of closest bed
open (PB, $proximityBed) || die "Can't open $proximityBed: $!\n";
while (<PB>) {
    $count++;
    chomp;
    my $line = $_;
    my @w = split(/\t/, $_);
    my $chr = $w[0];
    my $start = $w[1];
    my $stop = $w[2];
    my $entry = $chr . '_' . $start . '_' . $stop;
    next if ($entries{$entry});
    $entries{$entry} = 1;
    my $ampid = "OffRegion_" . $count;
    my $primerFwd = $count . "_fwd";
    my $primerRev = $count . "_rev";
    my $crossFwd = $alignment{$primerFwd}{hits}[0];
    my $crossRev = $alignment{$primerRev}{hits}[0];
    my $cigarFwd = $alignment{$primerFwd}{cigar};
    my $cigarRev = $alignment{$primerRev}{cigar};
    if (!$crossFwd || length($crossFwd) == 0) {
	$crossFwd = 'NA';
	$cigarFwd = 'NA';
    }
    if (!$crossRev || length($crossRev) == 0) {
	$crossRev = 'NA';
	$cigarRev = 'NA';
    }
    print OFFTARGETFLANKBEDALIGNED "$line\t$crossFwd\t$cigarFwd\t$crossRev\t$cigarRev\n";
    print OFFTARGETBEDNAMED "$chr\t$start\t$stop\t$ampid\t$crossFwd\t$cigarFwd\t$crossRev\t$cigarRev\n";
}
close (PB);
close (OFFTARGETFLANKBEDALIGNED);
close (OFFTARGETBEDNAMED);

$date = localtime;
print "Finished reading in proximity bed file in ion_offTargetCreateFlankBed.pl at $date\n";

$cmd = "mv $offTargetFlankBedAligned $proximityBed";
runAndCheckCommand(1);

$cmd = "cut -f1-4 $offTargetBedNamed > $offTargetBed";
runAndCheckCommand(1);

sub runAndCheckCommand {
    my $con = shift(@_);
    print($cmd, "\n");
    my $retVal += system($cmd);
    print "At ", `date`, "return value now is $retVal\n";
    if ($retVal != 0) {
	if ($con) {
	    $retVal = 0;
	    print "Program encountered problem after running $cmd, but keep going.\n";
	} else {
	    print "Program stop after running $cmd\n";
	    exit(0);
	}
    }
}
