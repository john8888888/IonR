#! /usr/bin/perl
use strict;
my $infile=$ARGV[0];
my $outfile=$ARGV[1];
open (OUT, ">$outfile") or die "Cannot open $outfile $!\n";
my $count=1;
my $track=0;
my $flag="--stranded=yes";
open(IN, $infile) or die "Cannot open $infile $!\n";
while(<IN>)
{
	chomp;
	my @array=split(/\t/,$_);
	next if ($_=~/^#/);
	if ($array[0]=~/track/)
	{
			$track++;
			next;
	}
#print STDERR "$array[5]\n";
	if ($array[5] eq '+' || $array[5] eq '-')
	{
		#$flag="--stranded=no";
		#$array[5]= '.';
	}
	else
	{
			$flag="--stranded=no";
	        $array[5]= '.';
	}
	if ($array[3] eq '')
	{
		$array[3]=$count;
	}
	if ($array[4] eq  '')
	{
			$array[4]=0;
	}
	print "$array[0]\tbed2gff\texon\t$array[1]\t$array[2]\t$array[4]\t$array[5]\t.\tgene_id \"$array[3]\"\n";
	$count++;
}
print OUT "$flag\n";
print STDERR "Number of tracks found=$track\n";
