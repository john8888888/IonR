#!/usr/bin/perl -w

use strict;
use warnings;

my $bamfile = $ARGV[0];
my $outputfile = $ARGV[1];
my $midfile = $ARGV[2];
#my $readLenFile = $ARGV[3];
my $readLenSum = $ARGV[3];
my $offTargetOn = $ARGV[4];
my $cutoffLength1 = $ARGV[5];
my $cutoffLength2 = $ARGV[6];
#my $cutoffLength3 = $ARGV[6] || die "Usage: original_bam output_file mid_sam_file_name read_length_summary 0 10 40\n";
my $cutoffLength3 = $ARGV[7] || die "Usage: original_bam output_file mid_sam_file_name read_length_summary 0 10 40 offTargetOn\n";

my $command = "samtools view $bamfile | cut -f1-10 > $midfile";
&runAndCheckCommand();
my $lineCt = 0;
my $ct0 = 0;
my $ct16 = 0;
my $ct4 = 0;
my $ctF0x4 = $lineCt - $ct4;
#20130415 new requested statistics, map quality score <4
my $cutoffMapQuality = 4;
my $ctLowerMQ = 0;
#count read length directly
my %lenDis;
my $cutoffCount1 = 0; #<=cutoffLength2                                                                           
my $cutoffCount2 = 0; #>cutoffLength2 && <cutoffLength3                                                          
my $cutoffCount3 = 0; #>cutoffLength3           

$offTargetOn = 1;
if ($offTargetOn) { #if not only essential analysis required
    open (MF, $midfile) || die "Can't open file: $!\n";
    while (<MF>) {
	chomp;
	$lineCt++;
	my @w = split(/\t/, $_);
	if ($w[1] == 0) {
	    $ct0++;
	} elsif ($w[1] == 16) {
	    $ct16++;
	} elsif ($w[1] == 4) {
	    $ct4++;
	} else {}


        #print RF "$w[0]\t", length($w[9]), "\n";
	my $l = length($w[9]);
	if ($l <= $cutoffLength2) {
	    $cutoffCount1++;
	} elsif ($l > $cutoffLength2 && $l <= $cutoffLength3) {
	    $cutoffCount2++;
	} else {
	    $cutoffCount3++;
	}
	

	$ctLowerMQ++ if ($w[4] < $cutoffMapQuality);
    }
    close (MF);
}

open (RF, ">$readLenSum") || die "Can't open file: $!\n";
#$lineCt++;
my $cutoffPerc1;# = $cutoffCount1 * 100 / $lineCt;
my $cutoffPerc2;# = $cutoffCount2 * 100 / $lineCt;
my $cutoffPerc3;# = $cutoffCount3 * 100 / $lineCt;
if ($offTargetOn) { #if not only essential analysis required
    $cutoffPerc1 = $cutoffCount1 * 100 / $lineCt;
    $cutoffPerc2 = $cutoffCount2 * 100 / $lineCt;
    $cutoffPerc3 = $cutoffCount3 * 100 / $lineCt;
} else {
    $cutoffPerc1 = -1;
    $cutoffPerc2 = -1;
    $cutoffPerc3 = -1;
}


print RF "$cutoffCount1\t";
printf RF "%.2f%%\t", $cutoffPerc1;
print RF 'length <= ', $cutoffLength1, "\n";
print RF "$cutoffCount2\t";
printf RF "%.2f%%\t", $cutoffPerc2;
print RF 'length > ', $cutoffLength1, ' and length <= ', $cutoffLength2, "\n";
print RF "$cutoffCount3\t";
printf RF "%.2f%%\t", $cutoffPerc3;
print RF 'length > ', $cutoffLength3, "\n";

close (RF);

#$lineCt = $lineCt - 3; #why 3 lines more?
$ctF0x4 = $lineCt - $ct4;
open (OF, ">$outputfile") || die "Can't open file: $!\n";
print OF "$lineCt\n$ct0\n$ct16\n$ct4\n$ctF0x4\n$ctLowerMQ\n";
close (OF);

sub runAndCheckCommand {
    my $con = shift(@_);
    print($command, "\n");
    my $retVal += system($command);
    print "At ", `date`, "return value now is $retVal\n";
    if ($retVal != 0) {
        if ($con) {
            $retVal = 0;
            print "Program encountered problem, but keep going.\n";
        } else {
            print "Program stop after running.\n";
	    exit(1);
        }
    }
}

