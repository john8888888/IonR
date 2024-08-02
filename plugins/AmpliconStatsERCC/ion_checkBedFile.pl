#!/usr/bin/perl -w

#check input bed file and remove duplicated ampliconID
#if duplicated, remove entry with POOL_ALIAS=Pool0; then pick the first entry

use strict;
use warnings;

my $inbed = $ARGV[0];
my $outbed = $ARGV[1] || die "Usage: input_bed_file output_bed_file\n";

my %entry;

open (IN, $inbed) || die "Can't open $inbed: $!\n";
while (<IN>) {
    next if (/^track/);

    chomp;
    my @f = split(/\s+/, $_);
    my $id = $f[3];

    if ($entry{$id}) {
	if ($entry{$id}{line} =~ /POOL_ALIAS=Pool0/i && !/POOL_ALIAS=Pool0/i) {
	    $entry{$id}{line} = $_;
	    $entry{$id}{chr} = $f[0];
	    $entry{$id}{start} = $f[1];
	    $entry{$id}{stop} = $f[2];
	} else {}
    } else {
	$entry{$id}{line} = $_;
	$entry{$id}{chr} = $f[0];
	$entry{$id}{start} = $f[1];
	$entry{$id}{stop} = $f[2];
    }
}
close (IN);

open (OUT, ">$outbed") || die "Can't open $outbed: $!\n";
foreach my $id (keys %entry) {
    print OUT "$entry{$id}{line}\n";
}
close (OUT);
