#!/usr/bin/perl -w

#to replace count5PrimeBase.sh

#$BASE_DIR/count5PrimeBase.sh $BASE_DIR $BAMFILE $REFERENCE >> ${OUTPUT_DIR}/$ONTARGETFIRST5BASEDAT;
#$BASE_DIR/count5PrimeBase.sh $BASE_DIR $OFFTARGETBAM $REFERENCE >> ${OUTPUT_DIR}/$OFFTARGETFIRST5BASEDAT;

use strict;
use warnings;

my $pluginDir = $ARGV[0];
my $bamFile = $ARGV[1];
my $reference = $ARGV[2];
my $out = $ARGV[3] || die "Usage: pluginDir bamFile reference outfile\n";

my $firstBaseBedFile = "first.5.prime.base.bed";
my $firstBaseBedFileFa = "$firstBaseBedFile.fa";

my $grandT = 0;
my $ttotal = 0;
my $atotal = 0;
my $ctotal = 0;
my $gtotal = 0;

my $direction;

my $cmd = "$pluginDir/ion_calcFirstBaseBed.py $bamFile f > f.$firstBaseBedFile";
runAndCheckCommand(1);
$cmd = "$pluginDir/ion_calcFirstBaseBed.py $bamFile F > F.$firstBaseBedFile";
runAndCheckCommand(1);
$cmd = "$pluginDir/fastaFromBed -s -name -tab -fi $reference -bed f.$firstBaseBedFile -fo stdout > f.$firstBaseBedFileFa";
runAndCheckCommand(1);
$cmd = "$pluginDir/fastaFromBed -s -name -tab -fi $reference -bed F.$firstBaseBedFile -fo stdout > F.$firstBaseBedFileFa";
runAndCheckCommand(1);

open (OUT, ">$out") || die "Can't open $out: $!\n";

my %ftotal;
open (FB, "f.$firstBaseBedFileFa") || die "Can't open f.$firstBaseBedFileFa: $!\n";
while (<FB>) {
    chomp;
    $ftotal{T}++ if (/\tT$/);
    $ftotal{A}++ if (/\tA$/);
    $ftotal{C}++ if (/\tC$/);
    $ftotal{G}++ if (/\tG$/);
}
close (FB);
my %Ftotal;
open (FB, "F.$firstBaseBedFileFa") || die "Can't open F.$firstBaseBedFileFa: $!\n";
while (<FB>) {
    chomp;
    $Ftotal{T}++ if (/\tT$/);
    $Ftotal{A}++ if (/\tA$/);
    $Ftotal{C}++ if (/\tC$/);
    $Ftotal{G}++ if (/\tG$/);
}
close (FB);

my $tt = $ftotal{T} + $Ftotal{T};
my $ta = $ftotal{A} + $Ftotal{A};
my $tc = $ftotal{C} + $Ftotal{C};
my $tg = $ftotal{G} + $Ftotal{G};
my $grand = $tt + $ta + $tc + $tg;

print OUT "Direction\tT\tA\tC\tG\n";
print OUT "Reverse\t$ftotal{T}\t$ftotal{A}\t$ftotal{C}\t$ftotal{G}\n";
print OUT "Forward\t$Ftotal{T}\t$Ftotal{A}\t$Ftotal{C}\t$Ftotal{G}\n";
print OUT "Total\t$tt\t$ta\t$tc\t$tg\n";
print OUT "Percent\t";
printf OUT "%.4f\t", $tt / $grand;
printf OUT "%.4f\t", $ta / $grand;
printf OUT "%.4f\t", $tc / $grand;
printf OUT "%.4f\n", $tg / $grand;
close (OUT);

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

