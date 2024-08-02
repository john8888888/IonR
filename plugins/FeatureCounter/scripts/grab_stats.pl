#! /usr/bin/perl
use strict;
my $infile=$ARGV[0]; # output from htseq-count
my $gff=$ARGV[1]; #input gff file
my $c_file=$ARGV[2]; # counts file
my $bam=$ARGV[3];
my $dir=$ARGV[4];
my $cmd="samtools flagstat $ARGV[3] > $dir/flagstats.txt";
#print "$cmd\n";
system($cmd);
open (TMP,"$dir/flagstats.txt") or die "Cannot open flagstats.txt in $dir $!\n";
my $mapped_reads;
while(<TMP>)
	{
		chomp;
		my $line=<TMP>;
		$line=<TMP>;
		my @map=split(/\s+/,$line);
		$mapped_reads=$map[0];
		last;
	}
my $count=0;
my ($no, $amb, $low, $not, $ali);
my %hash;
open (IN, $infile) or die "Cannot open $infile $!\n";
open (COU, ">$c_file") or die "Cannot open $c_file $!\n";
while(<IN>)
{
        chomp;
        my @array=split(/\t/,$_);

        if($_=~/^no_feature/)
        {
                $no=$array[1];
        }
        elsif ($_=~/^ambiguous/)
        {
                $amb=$array[1];
        }
        elsif ($_=~/^too_low_aQual/)
        {
                $low=$array[1];
        }
        elsif ($_=~/^not_aligned/)
        {
                $not=$array[1];
        }
        elsif ($_=~/^alignment_not_unique/)
        {
                $ali=$array[1];
        }
		else
		{
				$count+=$array[1];
				$hash{$array[0]}=$array[1];
		}

		
}
close(IN);
print COU "chr\tstart\tend\tfeature_id\tcounts\tRPKM\n";

open (GFF, $gff) or die "Cannot open $gff $!\n";
while(<GFF>)
{
        chomp;
        next if ($_=~/^#/);
	my @info=split(/\t/,$_);
        my $len=$info[4]-$info[3];
	$info[8]=~s/gene_id //;
	$info[8]=~s/\"//g;
	#print "$info[8]\n";
	my $RPKM; my $fRPKM;
	if ($count==0)
	{
		$RPKM=0;
	}
	else
	{
        	#print STDERR "$hash{$info[8]}\t$mapped_reads\n";
		$RPKM=(10**9*$hash{$info[8]})/($mapped_reads*$len);
		$fRPKM=rounded_to($RPKM,4);
	}
	print COU "$info[0]\t$info[3]\t$info[4]\t$info[8]\t$hash{$info[8]}\t$fRPKM\n";
}
my $total_reads=$mapped_reads;
#my $total_reads=$count+$no+$low+$ali+$amb+$not;
print "Total_reads: $total_reads\n";
print "Reads_on_features: $count\n";
print "No_feature: $no\n";
print "Ambiguous: $amb\n";
print "too_low_Quality: $low\n";
print "Not_Aligned: $not\n";
print "Alignment_not_unique: $ali\n";

sub rounded_to {    
    my ($round_me, $decimal_places) = @_;
    my $tmp = 10**$decimal_places;
    my $rounded = $round_me * $tmp;
    $rounded = int($rounded + 0.5 * ($rounded <=> 0));
    $rounded /= $tmp;
}
	
