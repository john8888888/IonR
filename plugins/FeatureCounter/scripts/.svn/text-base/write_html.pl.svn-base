#! /usr/bin/perl
use strict;

my $outputDirectory=$ARGV[0];
my $pluginpath=$ARGV[1];
my $urlRoot=$ARGV[2];
my $sampleName=$ARGV[3];

#parse output
my $total=`grep "Total_reads:" $outputDirectory/$sampleName/summary_stats.txt`;
my $count=`grep "Reads_on_features:" $outputDirectory/$sampleName/summary_stats.txt`;
my $no=`grep "No_feature:" $outputDirectory/$sampleName/summary_stats.txt`;
my $amb=`grep "Ambiguous:" $outputDirectory/$sampleName/summary_stats.txt`;
my $low=`grep "too_low_Quality:" $outputDirectory/$sampleName/summary_stats.txt`;
my $notAligned=`grep "Not_Aligned:" $outputDirectory/$sampleName/summary_stats.txt`;
my $notUnique=`grep "Alignment_not_unique:" $outputDirectory/$sampleName/summary_stats.txt`;

#print "$total\n";

chomp($total);
chomp($count);
chomp($no);
chomp($amb);
chomp($low);
chomp($notAligned);
chomp($notUnique);


$total=~s/Total_reads: //;
$count=~s/Reads_on_features: //;
$no=~s/No_feature: //;
$amb=~s/Ambiguous: //;
$low=~s/too_low_Quality: //;
$notAligned=~s/Not_Aligned: //;
$notUnique=~s/Alignment_not_unique: //;


print "Total is $total\n";


#if ($total ==0)
#{
#	$total=$no=$amb=$low=$notAligned=$notUnique=0
#}

#build html file from template
    my $command4 = "cat $pluginpath/template.html";
    my $htmlDocument = `$command4`;

    #print commands to txt file
    #open(OUT, ">$outputDirectory/$sampleName/commands.txt") || die "Couldn't write command text file\n";

    #print OUT "$command1\n";
    #print OUT "$command2\n";

    #close(OUT);

    #sub in our results
    $htmlDocument =~ s/#TOT_READS/$total/;
    $htmlDocument =~ s/#READS_OVL/$count/;
    $htmlDocument =~ s/#READS_NOT_OVL/$no/;
    $htmlDocument =~ s/#READS_AMB/$amb/;
    $htmlDocument =~ s/#READS_LOWQ/$low/;
    $htmlDocument =~ s/#READS_NOTA/$notAligned/;
    $htmlDocument =~ s/#READS_NOTU/$notUnique/;

    $htmlDocument =~ s/#URL_ROOT1/$urlRoot\/plugin_out\/FeatureCounter_out\/$sampleName\/feature_counts.txt/g;
    #$htmlDocument =~ s/#URL_ROOT2/$urlRoot\/plugin_out\/Assembler_out\/$sampleName\/$projectName\_assembly\/$projectName\_d_info\/$projectName/g;
    #$htmlDocument =~ s/#URL_ROOT3/$urlRoot\/plugin_out\/Assembler_out\/$sampleName/g;

    #print the html document
    my @nameToks = split(/\./, $sampleName);
    open(OUT, ">$outputDirectory/$nameToks[0]\.html") || die "Could not write sample html file\n";
    print OUT "$htmlDocument\n";
    close(OUT);

sub rounded_to {    
    my ($round_me, $decimal_places) = @_;
    my $tmp = 10**$decimal_places;
    my $rounded = $round_me * $tmp;
    $rounded = int($rounded + 0.5 * ($rounded <=> 0));
    $rounded /= $tmp;
}
