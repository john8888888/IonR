#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX;

my $inf1 = $ARGV[0];
my $inf2 = $ARGV[1];
my $out = $ARGV[2];
my $bias = $ARGV[3];
my $mergedBed = $ARGV[4];
my $offTargetOn = $ARGV[5];
my $chrStats = $ARGV[6] || die "Usage: for.depth.file rev.depth.file joined.depth.outfile bias.base.outfile merged.bed.file chr.stats.outfile offTargetOn\n";


open (INF1, $inf1) || die "Can't open $inf1: $!\n";
open (INF2, $inf2) || die "Can't open $inf2: $!\n";
open OUF, '>', $out or die "opening output file";
open (BIAS, ">$bias") || die "Can't open $bias: $!\n";
open (MBED, $mergedBed) || die "Can't open $mergedBed: $!\n";
open (CS, ">$chrStats") || die "Can't open $chrStats: $!\n";

my $biasedBase = 0;
my @in1;
my @in2;




my $numBaseInTarget = 0;
my $numBasesInTargetAutoChr = 0;
my %numBasesInTargetByChr;

while (<MBED>) {
    chomp;
    my @tokens = split("\t", $_);
    my $entryLen = $tokens[2] - $tokens[1];
    $numBaseInTarget += $entryLen;
    my $cid = uc($tokens[0]);
    if ($cid ne 'CHRX' && $cid ne 'CHRY') {
	$numBasesInTargetAutoChr += $entryLen;
    }
    $numBasesInTargetByChr{$cid} += $entryLen;
}

my $line; #the line to write into joined.depth.outfile
my $coverageLineCt = 0;
my %totalCoverageChr;
my %coverage;

# Prime the pump
&getrec1;
&getrec2;

my $preChr = $in1[0];
while (1) {
    $coverageLineCt++;
    print STDERR "$coverageLineCt processed at ", strftime("%m/%d/%Y %H:%M:%S\n", localtime), "\n" if ($coverageLineCt % 1000000 == 0);

    last if ($#in1<0 && $#in2<0);

    if ($#in1<0 || $#in2<0) {
	# Only one file is left...
	&write2 if $#in1<0;
	&write1 if $#in2<0;
    } else {
	my $write = decideWriteWhichOne($in1[0], $in1[1], $in2[0], $in2[1]);
	if ($write == 0) {
	    # Matching records, merge & write 'em
	    &writeboth;
	}
	elsif ($write == 1) {
	    # unmatched item in file 1, write it & get next rec
	    &write1;
	}
	else {
	    # unmatched item in file 2, write it & get next rec
	    &write2;
	}
    }

    my @tokens = split(/\t/, $line);
    my $cid = uc($tokens[0]);
    $totalCoverageChr{$cid} += $tokens[2];
    
    push @{ $coverage{$cid} }, $tokens[2];
}

print BIAS "$biasedBase bases with biased coverage (both strands >20x and one strand <3x)";
close (BIAS);

print STDERR "Preparing mean and 02mean  ", strftime("%m/%d/%Y %H:%M:%S\n", localtime), "\n";

exit if (!$offTargetOn);

my $totalCoverageAutoChr = 0;
my %coverageSummary;
foreach my $c (keys %totalCoverageChr) {
    my $mean = $totalCoverageChr{$c} / $numBasesInTargetByChr{$c};
    $coverageSummary{$c}{mean} = $mean;

    if ($c ne 'CHRX' && $c ne 'CHRY') {
	$totalCoverageAutoChr += $totalCoverageChr{$c};
    }
}
my $meanCoveragePerBaseAutoChr = $totalCoverageAutoChr / $numBasesInTargetAutoChr;
$coverageSummary{AutoChr}{mean} = $meanCoveragePerBaseAutoChr;

my %numBaseGT02MeanByChr;
my %percBaseGT02MeanByChr;
foreach my $c (keys %coverage) {
    foreach my $n (@{ $coverage{$c} }) {

	if ($n > $coverageSummary{$c}{mean} * 0.2) {
	    $coverageSummary{$c}{numBaseGT02Mean}++;

	    if ($c ne 'CHRX' && $c ne 'CHRY') {
		$coverageSummary{AutoChr}{numBaseGT02Mean}++;
	    }
	}
    }
}

foreach my $c (keys %coverageSummary) {
    if (uc($c) eq 'AUTOCHR') {
	$coverageSummary{$c}{percBaseGT02Mean} = $coverageSummary{$c}{numBaseGT02Mean} * 100 / $numBasesInTargetAutoChr;
    } else {
	$coverageSummary{$c}{percBaseGT02Mean} = $coverageSummary{$c}{numBaseGT02Mean} * 100 / $numBasesInTargetByChr{$c};
    }
}

print CS "chromosome\tmeanCoveragePerBaseByChr\tnumBaseGT02MeanByChr\tpercBaseGT02MeanByChr\n";
foreach my $c (sort sortChr (keys %coverage)) {
    printf CS "%s\t%2.2f\t%d\t%2.2f\n", ($c, $coverageSummary{$c}{mean}, $coverageSummary{$c}{numBaseGT02Mean}, $coverageSummary{$c}{percBaseGT02Mean});
}
printf CS "%s\t%2.2f\t%d\t%2.2f\n", ('AutoChr', $coverageSummary{AutoChr}{mean}, $coverageSummary{AutoChr}{numBaseGT02Mean}, $coverageSummary{AutoChr}{percBaseGT02Mean});
close(CS);

sub getrec1 {
    @in1 = ();
    if (!eof(INF1)) {
	(@in1) = split /\s+/, <INF1>;
	chomp $in1[0];
	chomp $in1[1];
	chomp $in1[2];
    }
}

sub getrec2 {
    @in2 = ();
    if (!eof(INF2)) {
	(@in2) = split /\s+/, <INF2>;
	chomp $in2[0];
	chomp $in2[1];
	chomp $in2[2];
    }
}

sub write1 {
    return if (!$in1[0]);
    $line = "$in1[0]\t$in1[1]\t$in1[2]\t$in1[2]\t0\n";
    print OUF $line; #"$in1[0]\t$in1[1]\t$in1[2]\t$in1[2]\t0\n";
    if ($in1[2] > 20) {
	$biasedBase++;
    }
    getrec1;
}

sub write2 {
    return if (!$in2[0]);
    $line = "$in2[0]\t$in2[1]\t$in2[2]\t0\t$in2[2]\n";
    print OUF $line; #"$in2[0]\t$in2[1]\t$in2[2]\t0\t$in2[2]\n";
    if ($in2[2] > 20) {
	$biasedBase++;
    }
    getrec2;
}

sub writeboth {
    my $sum = $in1[2] + $in2[2];
    $line = "$in1[0]\t$in1[1]\t$sum\t$in1[2]\t$in2[2]\n";
    print OUF $line; #"$in1[0]\t$in1[1]\t",$in1[2] + $in2[2],"\t$in1[2]\t$in2[2]\n";
    if (($in1[2] > 20 && $in2[2] < 3) || ($in2[2] > 20 && $in1[2] < 3)) {
	$biasedBase++;
    } 
    getrec1;
    getrec2;
}

sub decideWriteWhichOne {
    my ($c1, $p1, $c2, $p2) = @_;
    my $c1num = $c1;
    $c1num =~ s/chr//g;
    my $c2num = $c2;
    $c2num =~ s/chr//g;
    $c1num = 23 if ($c1num =~ /X/i);
    $c1num = 24 if ($c1num =~ /Y/i);
    $c1num = 25 if ($c1num =~ /M/i);
    $c2num = 23 if ($c2num =~ /X/i);
    $c2num = 24 if ($c2num =~ /Y/i);
    $c2num = 25 if ($c2num =~ /M/i);

    if ($c1num < $c2num) {
	return 1;
    } elsif ($c1num > $c2num) {
	return 2;
    } else {
	if ($p1 < $p2) {
	    return 1;
	} elsif ($p1 > $p2) {
	    return 2;
	} else {
	    return 0;
	}
    }
}

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
