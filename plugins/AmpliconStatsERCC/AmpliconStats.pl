#!/usr/bin/env perl
# Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved

use strict;
use warnings;
use FileHandle;
use Cwd;
use File::Temp;
use File::Basename;
use Getopt::Long;

my $opt = {
	"analysis-dir"  => undef,
	"out-dir"       => ".",
	"pre-out"       => undef,
	"pre-num"       => undef,
	"off-target"    => undef,
        "bed-file"      => undef,
	"bam-file"      => undef,
	"reference"     => undef,
	"web-bam"       => undef,
	"fastq"         => undef,
	"analysis-name" => undef,
	"library"       => undef,
	"run-name"      => undef,
	"base-name"     => "AmpliconStats",
	"version"       => undef,
	"help"          => 0,
	"just-html"     => 0,
};

GetOptions(
	"a|analysis-dir=s"  => \$opt->{"analysis-dir"}, 
	"o|out-dir=s"       => \$opt->{"out-dir"},
	"p|pre-out=s"       => \$opt->{"pre-out"},
	"u|pre-num=s"       => \$opt->{"pre-num"},
	"t|off-target=s"    => \$opt->{"off-target"},
        "b|bed-file=s"      => \$opt->{"bed-file"},
	"m|bam-file=s"      => \$opt->{"bam-file"},
        "f|reference=s"     => \$opt->{"reference"},
	"w|web-bam=s"       => \$opt->{"web-bam"},
	"q|fastq=s"         => \$opt->{"fastq"},
	"n|analysis-name=s" => \$opt->{"analysis-name"},
	"l|library=s"       => \$opt->{"library"},
	"r|run-name=s"      => \$opt->{"run-name"},
	"s|base-name=s"     => \$opt->{"base-name"},
	"v|version=s"       => \$opt->{"version"},
	"h|help"            => \$opt->{"help"},
	"j|just-html"       => \$opt->{"just-html"},
);

&usage() if(
	$opt->{"help"}
	|| !defined($opt->{"analysis-dir"})
	|| !defined($opt->{"bed-file"})
	|| !defined($opt->{"out-dir"})
	|| !defined($opt->{"bam-file"})
	|| !defined($opt->{"reference"})
);

sub usage () {
	print STDERR << "EOF";

	 usage: $0 [-h] --analysis-dir AnalysisDirPath --out-dir OutputDirPath --bed-file BedFilePath --bam-file BamFilePath --reference refereceFastaPath --web-bam BamFileULPath --fastq FastqFilePath
	     -a,--analysis-dir       : REQUIRED, directory with analysis results
	     -o,--out-dir            : REQUIRED, directory to write in
	     -p,--pre-out            : optional, should get previous results
	     -u,--pre-num            : optional, previous same name plugin output number
	     -t,--off-target         : optional, should turn on off-target analysis
	     -b,--bed-file           : REQUIRED, path to bed file
	     -m,--bam-file           : REQUIRED, path to bam file
	     -f,--reference          : REQUIRED, path to reference file
	     -w,--web-bam            : optional, path to url/web bam
	     -q,--fastq              : optional, path to fastq file
	     -n,--analysis-name      : optional, analysis name
	     -l,--library            : optional, library name
	     -r,--run-name           : optional, run name
	     -s,--base-name          : optional, base name, default AmpliconStats
	     -v,--version            : optional, version
	     -h,--help               : optional, print this help message
	     -j,--just-html          : optional, re-run to generate html page
EOF
	exit(0);
}


my $hostname = `hostname`;
print "HOSTNAME=$hostname\n";

#set variables
my $analysisDir = sprintf("%s",$opt->{"analysis-dir"});
my $plotSubDir = sprintf("%s",$opt->{"out-dir"});
my $targetBedFile = sprintf("%s",$opt->{"bed-file"});
my $bamOrig = sprintf("%s",$opt->{"bam-file"});
die "No bam file found. Exit." if (! -e $bamOrig);

my $webbam;
if ($opt->{"web-bam"}) {
    $webbam = sprintf("%s",$opt->{"web-bam"});
}
my $fastqOrig;
if ($opt->{"fastq"}) {
    $fastqOrig = sprintf("%s",$opt->{"fastq"});
}
my $genomeFasta = sprintf("%s",$opt->{"reference"});
my $version = $opt->{"version"};
if (!defined($version)){
    $version = "0.0.0";
}

my $analysisName = "myAnalysis";
my $libraryName = "hg19";
my $runName = "myRun";
if ($opt->{"analysis-name"}) {
    $analysisName = $opt->{"analysis-name"};
}
if ($opt->{"library"}) {
    $libraryName = $opt->{"library"};
}
if ($opt->{"run-name"}) {
    $runName = $opt->{"run-name"};
}
$analysisName .= '.AmpStats';

my $command = '';
my $retVal = 0;



#mkdir $opt->{"out-dir"} || die "$0: unable to make directory ".$opt->{"out-dir"}.": $!\n";
mkdir $plotSubDir || die "$0: unable to make directory $plotSubDir $!\n";
my ($tmpfile,$plugindir) = fileparse($0);

my $justhtml = 0;
if ($opt->{"pre-out"} eq 'Resume') {
    #find previous thisPlugin_out.xxxxxx directory
    my $preOutDir = &findPreviousOutputFolder($plotSubDir, $opt->{"pre-num"});
    print "preOutDir is $preOutDir\n";
    my $justhtmlFile = $preOutDir . "/return.value.txt";
    if (-e $justhtmlFile) {
	$command = "cp $preOutDir/* $plotSubDir/ 1>/dev/null 2>/dev/null";
	&runAndCheckCommand(1);
	$command = "cp $preOutDir/splitBeds/* $plotSubDir/ 1>/dev/null 2>/dev/null";
	&runAndCheckCommand(1);
	$justhtml = 1;
	print "Generating HTML only.\n";
    }
} else {
    $justhtml = 0;
}
my $offTargetOn = 0;
if ($opt->{"off-target"} eq 'On') {
    $offTargetOn = 1;
}

#exit(0);

for my $f ('SummaryTable.txt', 'results.json', 'AmpliconStats.html', 'summary_table.txt') {
    if (-e "$plotSubDir/$f") {
	$command = "mv $plotSubDir/$f $plotSubDir/$f.`date +%s`";
	&runAndCheckCommand(1);
    }
}

$command = "rm " . "$plotSubDir/*" . $analysisName . "*" . " 1>/dev/null 2>/dev/null";
&runAndCheckCommand(1);
#$command = "rm " . "$plotSubDir/*" . 'thumbnail.png';
#&runAndCheckCommand(1);
$command = "rm $plotSubDir". "/return.value.txt" . " 1>/dev/null 2>/dev/null";
&runAndCheckCommand(1);
$command = "rm $plotSubDir". '/*first.5.prime.base*' . " 1>/dev/null 2>/dev/null";
&runAndCheckCommand(1);

#used for total base number calculation, removing overlapping ones
my @temp = split(/\//, $targetBedFile);
my $targetBedFileCopy = $targetBedFile;
my $mergedBed = $plotSubDir . '/' . $analysisName . '.' . $temp[-1] . '.merged'; 
my $copyBed = $plotSubDir . '/' . $analysisName . '.' . $temp[-1] . '.copy'; 
$command = "cp $targetBedFile $copyBed 1>/dev/null 2>/dev/null";
&runAndCheckCommand();
$command = "$plugindir/mergeBed -i $targetBedFile > $mergedBed";
&runAndCheckCommand();
my $ampliconFasta = $plotSubDir . '/../' . $temp[-1] . '.fasta';

#remove duplicate amplicon ID
my $checkedBed = $plotSubDir . '/' . $analysisName . '.' . $temp[-1] . '.checked'; 
$command = "perl $plugindir/ion_checkBedFile.pl $targetBedFile $checkedBed";
if ($offTargetOn) {
    &runAndCheckCommand();
    $targetBedFile = $checkedBed;
} else {
    #copied to here
}

#$command = "cp $targetBedFile $copyBed";
#&runAndCheckCommand();

my $samtools = $plugindir."/samtools";

#generate directional depth
my $hg19ontargetCoverage = $analysisName . ".hg19_ontarget_coverage.txt";
#my $hg19ontargetCoverageShort = "hg19_ontarget_coverage.txt.short";
#my $hg19ontargetCoverageStrandF = "hg19_ontarget_coverage.strand.F.txt";
#my $hg19ontargetCoverageStrandR = "hg19_ontarget_coverage.strand.R.txt";
my $readLenName = "$plotSubDir/$analysisName" . ".readLen.txt";
my $readLenSum = "$plotSubDir/$analysisName" . ".readLen.summary.txt";
my $posErrorName = "$plotSubDir/$analysisName" . ".position_error_summary.txt";
my $ampsum = "$plotSubDir/$analysisName" . "_ampliconSummary.csv";
my $midSamFile = "$plotSubDir/$analysisName" . ".mid.sam";
my $midF0x4SamFile = "$plotSubDir/$analysisName" . ".mid.F0x4.sam";
my $readInfoFile = "$plotSubDir/$analysisName" . ".readInfo.txt";
my $biasedBaseFile = "$plotSubDir/$analysisName" . ".biasedBase.txt";

#if (! $opt->{"just-html"}) {
if (! $justhtml) {
    chdir($opt->{"out-dir"});
    #&genReadSummary($bamOrig, "$plotSubDir/$analysisName" . ".readInfo.txt", $midSamFile);
    #$command = "samtools view $bamOrig | awk '{print \$1\"\t\"length(\$10)}' > $readLenName";
    #above two lines were replaced by ion_genReadSummary.pl
    #$command = "perl $plugindir/ion_genReadSummary.pl $bamOrig $readInfoFile $midSamFile $readLenName $readLenSum 0 10 40";
    #$command = "perl $plugindir/ion_genReadSummary.pl $bamOrig $readInfoFile $midSamFile $readLenSum 0 10 40";
    $command = "perl $plugindir/ion_genReadSummary.pl $bamOrig $readInfoFile $midSamFile $readLenSum $offTargetOn 0 10 40";
    &runAndCheckCommand();
    #&getReadLengthSummary($readLenName, $readLenSum, 0, 10, 40);

    if(-e $targetBedFile){
	# generate the amplicon summary table
	#my $genomeFasta = "/results/referenceLibrary/tmap-f2/$libraryName/$libraryName.fasta";
	if (-e $ampliconFasta) {
	    #no need run this for every barcode if the same bed file is used
	} else {
	    $command = "$plugindir/fastaFromBed -name -fi $genomeFasta -bed $targetBedFile -fo $ampliconFasta";
	    &runAndCheckCommand();
	}

	#run ion_parseBAM.sh
	#$command = "$plugindir/ion_parseBAM.sh $bamOrig $targetBedFileCopy $libraryName $analysisName $ampliconFasta $plugindir $plotSubDir $hg19ontargetCoverage $ampsum $readInfoFile $midF0x4SamFile";
	my $chrstats = "$plotSubDir/$analysisName.base.coverage.by.chr.tab";
	#$command = "$plugindir/ion_parseBAM.sh $bamOrig $targetBedFileCopy $libraryName $analysisName $ampliconFasta $plugindir $plotSubDir $hg19ontargetCoverage $ampsum $readInfoFile $midF0x4SamFile $biasedBaseFile $mergedBed $chrstats";
	$command = "$plugindir/ion_parseBAM.sh $bamOrig $targetBedFileCopy $libraryName $analysisName $ampliconFasta $plugindir $plotSubDir $hg19ontargetCoverage $ampsum $readInfoFile $midF0x4SamFile $biasedBaseFile $mergedBed $chrstats $offTargetOn";
	#$command = "$plugindir/ion_parseBAM.sh $bamOrig $targetBedFileCopy $libraryName $analysisName $ampliconFasta $plugindir $plotSubDir $hg19ontargetCoverage $ampsum $readInfoFile $midF0x4SamFile $biasedBaseFile";
	&runAndCheckCommand();

	if (-e "$plotSubDir/$hg19ontargetCoverage"){
	    #$command = "python $plugindir/ion_chromosomeReporter.py $bamOrig $targetBedFile $plotSubDir/$hg19ontargetCoverage > $plotSubDir/$analysisName"."_chromosomeSummary.csv";
	    if ($offTargetOn) {
	    my $chrReportFile = "$plotSubDir/$analysisName"."_chromosomeSummary.csv";
	    $command = "perl $plugindir/ion_position_error_and_chromosome_report.pl $targetBedFile $bamOrig $plotSubDir/$hg19ontargetCoverage $posErrorName $chrReportFile";
	    &runAndCheckCommand();
	    }
	    #generate positional error data
	    #$command = "perl $plugindir/ion_position_error.pl $targetBedFile $bamOrig $samtools > $posErrorName";
	    #&runAndCheckCommand();
	}
    } else {
	my $msg = 'No specified target bed file could be found.';
	&writeUnsuitableHtml($opt,$analysisName, "AmpliconStats.html", $plotSubDir, "Amplicon General Analysis", $version, $msg);
	exit(0);
    }

    if ($offTargetOn) {
	$command = "$plugindir/offTarget.sh $bamOrig $targetBedFile $plotSubDir $plugindir $genomeFasta $analysisName $ampsum";
	&runAndCheckCommand();
    }

    if (-e "$plotSubDir/$analysisName.targetStats.txt"){
	my $ampSumFile = $plotSubDir . '/' . $analysisName . '_ampliconSummary.csv';
	#hsm loci priority data
	my $hsmPri1File;
	my $hsmPri9File;

	if ($targetBedFile =~ /400_hsm_v12_1_seq/) {
	    $hsmPri1File = $plotSubDir . "/hsm_stats_norm_pri1.txt";
	    $hsmPri9File = $plotSubDir . "/hsm_stats_norm_pri9.txt";
	    my $hotspotBed = "$plugindir/bedfiles/HotSpots_1.0_Ion_AmpliSeq_Cancer.bed";
	    my $priorityBed = "$plugindir/bedfiles/400_hsm_v12_1.extended.bed";
	    my $priCut = 500;

	    &getPriorityData($bamOrig, $hotspotBed, $priorityBed, $targetBedFile, $priCut, $hsmPri1File, $hsmPri9File);

	} else {
	    $hsmPri1File = $plotSubDir . "/doesnotexist.txt";
	    $hsmPri9File = $plotSubDir . "/doesnotexist.txt";
	}

	#move calculation to ion_join.coverage.depth.pl in ion_parseBAM.sh
	#$command = "python $plugindir/ion_getChromosomalBaseStats.py $plotSubDir/$hg19ontargetCoverage $mergedBed > $plotSubDir/$analysisName.base.coverage.by.chr.tab";
	#&runAndCheckCommand();

	#$command = "python $plugindir/ion_generateStats.py $plotSubDir/$analysisName.targetStats.txt $targetBedFile $plotSubDir/$hg19ontargetCoverage $ampSumFile $mergedBed $readLenSum $analysisDir/bfmask.stats $posErrorName $plotSubDir $plotSubDir/results.json $plotSubDir/summary_table.txt $hsmPri1File $hsmPri9File $plotSubDir/$analysisName.base.coverage.by.chr.tab $biasedBaseFile $readInfoFile";
	$command = "python $plugindir/ion_generateStats.py $plotSubDir/$analysisName.targetStats.txt $targetBedFile $plotSubDir/$hg19ontargetCoverage $ampSumFile $mergedBed $readLenSum $posErrorName $plotSubDir $plotSubDir/results.json $plotSubDir/summary_table.txt $hsmPri1File $hsmPri9File $plotSubDir/$analysisName.base.coverage.by.chr.tab $biasedBaseFile $readInfoFile";
	&runAndCheckCommand();
	#generate thumbnail image for all png files
	$command = "python $plugindir/png.to.thumbnail.py";
	&runAndCheckCommand(1);
    }
}# else {
my $jhTemp = $justhtml;
$justhtml = 0; #reset
    $command = "mkdir $plotSubDir/splitBeds/";
    &runAndCheckCommand(1);
    $command = "mv $plotSubDir/split*bed $plotSubDir/splitBeds/";
    &runAndCheckCommand(1);
$justhtml = $jhTemp; #reset
#}

# Cleanup

# Write the html output
#    my($opt, $retVal, $analysisName, $htmlFile, $htmlDir, $plotDir, $customTitle, $version, $bamOrig, $hg19ontargetCoverage) = @_;
&writeHtml($opt, $retVal, $analysisName, "AmpliconStats.html", $plotSubDir, $plotSubDir, "Amplicon Statistics", $version, $bamOrig, $hg19ontargetCoverage);

my $retValFile = "$plotSubDir/return.value.txt";
open (RV, ">$retValFile") || die "Can't open $retValFile: $!\n";
print RV "$retVal\n";
close (RV) || die "Can't close $retValFile: $!\n";
#exit($retVal);
exit(0);

sub writeTabFile {
    my($htmlFh, $filename, $sortable) = @_;

    my $tableIn = FileHandle->new($filename) || die "$0: Can't open file. $!\n";
    my @pngPrefixes;
    my $count = 0;
    if (!$sortable) {
	print $htmlFh "<table border=1 cellpadding=6>";
    } else {
	print $htmlFh "<table  id='metrics_tab' name='metrics_tab' class='sortable' border=1 cellpadding=6 width=100%>";
    }

    while (my $line = <$tableIn>) { 
	chomp($line);
	my @entries = split /\t/,$line;
	print $htmlFh "<tr>";
	for (my $i = 0; $i < scalar(@entries); $i++) {
	    print $htmlFh "<td>$entries[$i]</td>";
	}
	print $htmlFh "</tr>\n";
	$pngPrefixes[$count++] = $entries[0];
    }
    print $htmlFh "</table>";

}

sub writeCSVFile {
    my($htmlFh, $filename) = @_;

    my $tableIn = FileHandle->new($filename) || die "$0: Can't open file. $!\n";
    my $header = <$tableIn>; 
    chomp($header);
    my @colHeader = split /,/,$header;
    my $count = 0;
    print $htmlFh "<table border=1 cellpadding=6>";
    print $htmlFh "<tr>";
    for (my $i = 0; $i < scalar(@colHeader); $i++) {
	print $htmlFh "<th>$colHeader[$i]</th>";
    }
    print $htmlFh "</tr>\n";
    while (my $line = <$tableIn>) { 
	chomp($line);
	my @entries = split /,/,$line;
	print $htmlFh "<tr>";
	for (my $i = 0; $i < scalar(@entries); $i++) {
	    print $htmlFh "<td>$entries[$i]</td>";
	}
	print $htmlFh "</tr>\n";
    }
    print $htmlFh "</table>";

}

sub writeTextFile {
    my($htmlFh, $filename) = @_;

    my $tableIn = FileHandle->new($filename) || die "$0: Can't open file. $!\n";
    while (my $line = <$tableIn>) { 
	chomp($line);
	print $htmlFh "$line<br>";
    }
    print $htmlFh "<br><br>";
}


sub genReadSummary {
   my($bamfile, $outputfile, $midfile) = @_;

    $command = "$samtools view $bamfile | cut -f1-10 > $midfile";
    &runAndCheckCommand();
    my $lineCt = 0;
    my $ct0 = 0;
    my $ct16 = 0;
    my $ct4 = 0;
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
    }
    close (MF);
    open (OF, ">$outputfile") || die "Can't open file: $!\n";
    print OF "$lineCt\n$ct0\n$ct16\n$ct4\n";
    close (OF);

    #$retVal += system("$samtools view $bamfile | wc -l > $outputfile");
    #$retVal += system("$samtools view $bamfile | cut -f2 | grep 0 | wc -l >> $outputfile");
    #$retVal += system("$samtools view $bamfile | cut -f2 | grep 16 | wc -l >> $outputfile");
    #$retVal += system("$samtools view $bamfile | cut -f2 | grep 4 | wc -l >> $outputfile");
    $retVal += system("$samtools view -F 0x4 $bamfile | wc -l >> $outputfile");
}

=hold
sub getReadLengthSummary {
    my ($readLengthFile, $readLengthOutput, $cutoffLength1, $cutoffLength2, $cutoffLength3) = @_;
    my %readLength;
    my $cutoffCount1 = 0; #<=cutoffLength2
    my $cutoffCount2 = 0; #>cutoffLength2 && <cutoffLength3
    my $cutoffCount3 = 0; #>cutoffLength3
    my $totalRead = 0;

    open (RLOUT, ">$readLengthOutput") || die "Can't open $readLengthOutput: $!\n";
    open (RL, $readLengthFile) || die "Can't open $readLengthFile: $!\n";
    <RL>;
    while (<RL>) {
        chomp;
        my ($id, $len) = split(/\t/, $_);
        $readLength{$len}++;
        $totalRead++;
    }
    close (RL) || die "Can't close $readLengthFile: $!\n";

    for my $l (sort {$a <=> $b} keys %readLength) {
        $cutoffCount1 += $readLength{$l} if ($l <= $cutoffLength2);
        $cutoffCount2 += $readLength{$l} if ($l > $cutoffLength2 && $l <= $cutoffLength3);
        $cutoffCount3 += $readLength{$l} if ($l > $cutoffLength3);
    }

    print RLOUT "$cutoffCount1\t";
    printf RLOUT "%4.2f", $cutoffCount1/$totalRead*100;
    print RLOUT "%\tlength <= $cutoffLength2\n";
    print RLOUT "$cutoffCount2\t";
    printf RLOUT "%4.2f", $cutoffCount2/$totalRead*100;
    print RLOUT "%\tlength > $cutoffLength2 and length <= $cutoffLength3\n";
    print RLOUT "$cutoffCount3\t";
    printf RLOUT "%4.2f", $cutoffCount3/$totalRead*100;
    print RLOUT "%\tlength > $cutoffLength3\n";

    print RLOUT "$totalRead\tTotal reads\n";
    close (RLOUT) || die "Can't close $readLengthOutput: $!\n";

}
=cut

sub writeUnsuitableHtml {
    my($opt,$analysisName, $htmlFile, $htmlDir, $customTitle, $version, $content) = @_;
    my $cwd = &Cwd::getcwd(); #save so we can go back
    chdir $htmlDir; 
    #this html page at top level

    my $htmlFh = FileHandle->new("> $htmlFile") || die "$0: problem writing $htmlFile: $!\n";
    my $plotFile = ""; 

    my $title = $opt->{"base-name"};

    &writeHtmlHeader($title, $customTitle."<br>Version ".$version."<br>".$analysisName,$htmlFh);
    print $htmlFh "<h1>", $content, "</h1>";
    finishHtml($htmlFh);
    chdir $cwd; #go back
    return;
}

#&writeHtml($opt, $retVal, $analysisName, "AmpliconStats.html", $plotSubDir, $plotSubDir, "Amplicon Statistics", $version, $bamOrig, $hg19ontargetCoverage);
sub writeHtml {
    my($opt, $retVal, $analysisName, $htmlFile, $htmlDir, $plotDir, $customTitle, $version, $bamOrig, $hg19ontargetCoverage) = @_;
    my $cwd = &Cwd::getcwd(); #save so we can go back
    chdir $htmlDir; 
    #this html page at top level
    $command = "ln -sf $plugindir/js .";
    runAndCheckCommand(1);
    $command = "ln -sf $plugindir/css .";
    runAndCheckCommand(1);
    $command = "ln -sf $plugindir/export .";
    runAndCheckCommand(1);

    my $htmlFh = FileHandle->new("> $htmlFile") || die "$0: problem writing $htmlFile: $!\n";
    my $plotFile = ""; 

    my $title = $opt->{"base-name"};

    #cp bedfile etc for download
    my @temp = split(/\//, $targetBedFileCopy);
    my $targetBedFileName = $temp[-1];
    $command = "cp -f $targetBedFileCopy ./";
    &runAndCheckCommand(1);


    &writeHtmlHeader($title, $customTitle."<br>Version ".$version."<br>".$analysisName, $htmlFh);
    if ($retVal != 0) { #should be already caught in previous steps
	print $htmlFh "Encountered error while running the plugin.\n";
	finishHtml($htmlFh);
	close($htmlFh);
 	#exit($retVal);
 	exit(0);
    }
    print $htmlFh "<center>Reference: $libraryName<br>Target file <a href=\"$targetBedFileName\">$targetBedFileName</a><br><br></center>\n";

    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    $mon++;
    $year =~ s/11/1/; #hard code here to fix year number 113 to 2013
    print $htmlFh "<center>Running date: $mon/$mday/$year<br><br></center>\n";


    #read ampliconSummary.csv and summary_table.txt to get values of top table
    my $ampSum = $analysisName."_ampliconSummary.csv";
    my $sumTab = "summary_table.txt";
    my $sumCov = "coverage_table.txt";
    if ((-e $sumTab) && (-e $ampSum)){
	print $htmlFh "<center>\n";
	print $htmlFh "<table border=1 cellpadding=6>\n";

	my $field;
	$field = 'Number of targets';
	my $numOfTarget = `grep "$field" $sumTab | cut -f2`;
	$field = 'Per base accuracy';
	my $perBaseAcc = `grep "$field" $sumTab | cut -f2`;
	$field = 'Percent greater than 0.2 mean reads per base';
	my $perc02Mean = `grep "$field" $sumTab | cut -f2`;
	$field = 'Percent no strand bias of all bases';
	my $percBaseNoBias = `grep "$field" $sumTab | cut -f2`;
	#$field = 'Percent all reads on target';
	$field = 'Percent mapped reads on target';
	my $percReadsOnTarget = `grep "$field" $sumTab | cut -f2`;
	$field = 'Target coverage at 20x - norm 100';
	my $percTargetCoverage20xNorm100 = `grep "$field" $sumTab | cut -f2`;
	#my $percTargetCoverage20xNorm100 = `grep -i normalized $sumCov | cut -f4`;
	my $passNumOfTarget = `cut -f12 -d, $ampSum | grep -v Pass | grep -c 1`;
	my $passNumOfTarget001 = `cut -f28 -d, $ampSum | grep -v Pass | grep -c 1`;
	#chomp($numOfTarget, $perBaseAcc, $perc02Mean, $percBaseNoBias, $percReadsOnTarget, $passNumOfTarget, $percTargetCoverage20xNorm100);
	chomp($numOfTarget, $perBaseAcc, $perc02Mean, $percBaseNoBias, $percReadsOnTarget, $passNumOfTarget, $passNumOfTarget001, $percTargetCoverage20xNorm100);
	my $assayConvRate = sprintf("%.2f%%", ($passNumOfTarget/$numOfTarget)*100);
	my $assayConvRate001 = sprintf("%.2f%%", ($passNumOfTarget001/$numOfTarget)*100);

        print $htmlFh "<tr><th align='left'>Number of targets</th><td>$numOfTarget</td></tr>\n";
        print $htmlFh "<tr><th align='left'>Per Base Accuracy</th><td>$perBaseAcc</td></tr>\n";
        print $htmlFh "<tr><th align='left'>Percent greater than 0.2 mean reads per base</th><td>$perc02Mean</td></tr>\n";
        print $htmlFh "<tr><th align='left'>Base without Strand Bias</th><td>$percBaseNoBias</td></tr>\n";
        #print $htmlFh "<tr><th align='left'>Percent of reads on target</th><td>$percReadsOnTarget</td></tr>\n";
        print $htmlFh "<tr><th align='left'>Percent of mapped reads on target</th><td>$percReadsOnTarget</td></tr>\n";
        #print $htmlFh "<tr><th align='left'>Assay Conversion Rate</th><td>$assayConvRate</td></tr>\n";
        print $htmlFh "<tr><th align='left'>Assay Conversion Rate at 0.2*mean</th><td>$assayConvRate</td></tr>\n";
        print $htmlFh "<tr><th align='left'>Assay Conversion Rate at 0.01*mean</th><td>$assayConvRate001</td></tr>\n";
        print $htmlFh "<tr><th align='left'>Target coverage at 20x - norm 100</th><td>$percTargetCoverage20xNorm100</td></tr>\n";

	print $htmlFh "</table>\n";
	print $htmlFh "</center>\n";
    }

    print $htmlFh "<br><br>\n";

    print $htmlFh "<br><h2>Table of Contents</h2>\n";
    print $htmlFh "<br><a href=\"\#DOWNLOADS\">Downloads</a>\n";
    print $htmlFh "<br><a href=\"\#SUMMARY\">Summary Table</a>\n";
    print $htmlFh "<br><a href=\"\#NumReads\">Amplicon Representation</a>\n";
    print $htmlFh "<br><a href=\"\#StatPlots\">Amplicon Statistics</a>\n";
    #print $htmlFh "<br><a href=\"\#ErrPos\">Error By Position</a>\n";
    print $htmlFh "<br><a href=\"\#OffTarget\">Off Target Analysis</a>\n";
    print $htmlFh "<br><a href=\"\#ImmediateBase\">Immediate 5-prime Base Analysis</a>\n";


    #download section
    print $htmlFh "<a name=\"DOWNLOADS\"><h3>Downloads</h3></a>\n";
    my @work = split(/\//, $bamOrig);
    my $bamOrigName = $work[-1];
    if (-e $bamOrigName) {
	print $htmlFh "      <a href=\"$bamOrigName\">HG19 Bam</a> (link to the bam file mapping all reads to HG19)<br>\n";
	print $htmlFh "      <a href=\"$bamOrigName\.bai\">HG19 Bam Index</a> (link to the bam index file mapping all reads to HG19)<br>\n";
    } elsif (-e "../../$bamOrigName") {
	print $htmlFh "      <a href=\"../../$bamOrigName\">HG19 Bam</a> (link to the bam file mapping all reads to HG19)<br>\n";
	print $htmlFh "      <a href=\"../../$bamOrigName\.bai\">HG19 Bam Index</a> (link to the bam index file mapping all reads to HG19)<br>\n";
    } else {
	my $webbamtemp = $webbam;
	$webbamtemp =~ s/results\d*\/analysis\///;
	print $htmlFh "      <a href=\"$webbamtemp\">HG19 Bam</a> (link to the bam file mapping all reads to HG19)<br>\n";
	print $htmlFh "      <a href=\"$webbamtemp\.bai\">HG19 Bam Index</a> (link to the bam index file mapping all reads to HG19)<br>\n";
    }

    $plotFile = sprintf("%s_ampliconSummary.csv",$analysisName);
    print $htmlFh "      <a href=\"$plotFile\">Amplicon Summary</a><br>\n" if(-e $plotFile);

    @temp = split(/\//, $targetBedFile);
    $targetBedFileName = $temp[-1];
    print $htmlFh "      <a href=\"$targetBedFileName\">Amplicon bed file ($targetBedFileName)</a><br>\n" if(-e $targetBedFileName);

    $plotFile = sprintf("%s", $hg19ontargetCoverage);
    print $htmlFh "<a href=\"$plotFile\">On Target Coverage Data</a>  (generated by samtools depth -b) <br>\n" if(-e $plotFile);

    #summary table
    print $htmlFh "<a name=\"SUMMARY\"><h3>Summary Table</h3></a>\n";
    print $htmlFh "<center>\n";
    $plotFile = sprintf("summary_table.txt");
    print $htmlFh "<a href=\"$plotFile\">Download</a><br>\n" if(-e $plotFile);
    &writeTabFile($htmlFh, $plotFile) if (-e $plotFile);
    print $htmlFh "</center>\n";

=take coverage table off
   #coverage table
    print $htmlFh "<a name=\"COVERAGE\"><h3>Coverage Table</h3></a>\n";
    print $htmlFh "<center>\n";
    $plotFile = sprintf("coverage_table.txt");
    print $htmlFh "<a href=\"$plotFile\">Download</a><br>\n" if(-e $plotFile);
    &writeTabFile($htmlFh, $plotFile) if (-e $plotFile);
    print $htmlFh "</center>\n";
=cut

#$analysisName.base.coverage.by.chr.ta
    # the link to base coverage by chromosome summary table.
    $plotFile = sprintf("%s.base.coverage.by.chr.tab",$analysisName);
    print $htmlFh "      <h3>Base Coverage Summary by Sex/Auto-Chromosome</h3>\n" if(-e $plotFile);
    print $htmlFh "<center>\n";
    print $htmlFh "      <a href=\"$plotFile\">Download</a><br>\n" if(-e $plotFile);
    &writeTabFile($htmlFh, $plotFile) if (-e $plotFile);
    print $htmlFh "</center><br>\n";

    # the link to the chromosome summary table.
    $plotFile = sprintf("%s_chromosomeSummary.csv",$analysisName);
    print $htmlFh "      <h3>Chromosome Coverage Summary</h3>\n" if(-e $plotFile);
    print $htmlFh "<center>\n";
    print $htmlFh "      <a href=\"$plotFile\">Download</a><br>\n" if(-e $plotFile);
    &writeCSVFile($htmlFh, $plotFile) if (-e $plotFile);
    print $htmlFh "</center><br>\n";

    # Plot Amplicon Represenation
    print $htmlFh "<a name=\"NumReads\"><h3>Amplicon Representation</h3></a>\n";
    print $htmlFh "<table border=1 cellpadding=3 width=100%>\n";
    print $htmlFh "    <tr>\n";
    $plotFile = sprintf("%s_representation.png", $analysisName);
    print $htmlFh "<td><a href=\"$plotFile\"><img width=100% src=\"$plotFile\"/></a></td>\n" if(-e $plotFile);
    $plotFile = sprintf("%s_log_representation.png", $analysisName);
    print $htmlFh "<td><a href=\"$plotFile\"><img width=100% src=\"$plotFile\"/></a></td>\n" if(-e $plotFile);
    $plotFile = sprintf("%s_repBias.png", $analysisName);
    print $htmlFh "<td><a href=\"$plotFile\"><img width=100% src=\"$plotFile\"/></a></td>\n" if(-e $plotFile);
    print $htmlFh "    </tr>\n";

    $plotFile = sprintf("%s_sucFailvsGC.png", $analysisName);
    print $htmlFh "<td><a href=\"$plotFile\"><img width=100% src=\"$plotFile\"/></a></td>\n" if(-e $plotFile);
    $plotFile = sprintf("%s_sucFailRatevsGC.png", $analysisName);
    print $htmlFh "<td><a href=\"$plotFile\"><img width=100% src=\"$plotFile\"/></a></td>\n" if(-e $plotFile);
    $plotFile = sprintf("predicted_coverage.png");
    print $htmlFh "<td><a href=\"$plotFile\"><img width=100% src=\"$plotFile\"/></a></td>\n" if(-e $plotFile);
    print $htmlFh "    </tr>\n";
    print $htmlFh "</table>\n";

    # Plot Amplicon Statistics
    print $htmlFh "<a name=\"StatPlots\"><h3>Amplicon Statistics</h3></a>\n";
    print $htmlFh "<table border=1 cellpadding=3 width=100%>\n";
    print $htmlFh "    <tr>\n";
    $plotFile = sprintf("%s_gc_fwd.png", $analysisName);
    print $htmlFh "<td><a href=\"$plotFile\"><img width=100% src=\"$plotFile\"/></a></td>\n" if(-e $plotFile);
    $plotFile = sprintf("%s_gc_rev.png", $analysisName);
    print $htmlFh "<td><a href=\"$plotFile\"><img width=100% src=\"$plotFile\"/></a></td>\n" if(-e $plotFile);
    $plotFile = sprintf("%s_len_fwd.png", $analysisName);
    print $htmlFh "<td><a href=\"$plotFile\"><img width=100% src=\"$plotFile\"/></a></td>\n" if(-e $plotFile);
    $plotFile = sprintf("%s_len_rev.png", $analysisName);
    print $htmlFh "<td><a href=\"$plotFile\"><img width=100% src=\"$plotFile\"/></a></td>\n" if(-e $plotFile);
    print $htmlFh "    </tr>\n";
    print $htmlFh "</table>\n";

    # Plot Off Target GC Boxplot
    print $htmlFh "<a name=\"GCPlots\"><h3>GC Comparison</h3></a>\n";
    print $htmlFh "<table border=1 cellpadding=3 width=25%>\n";
    print $htmlFh "    <tr>\n";
    $plotFile = sprintf("%s_gc.on.vs.off.target.png", $analysisName);
    print $htmlFh "<td><a href=\"$plotFile\"><img width=100% src=\"$plotFile\"/></a></td>\n" if(-e $plotFile);
    $plotFile = sprintf("%s_gc.on.vs.off.target.png_nonexistant", $analysisName);
    print $htmlFh "    </tr>\n";
    print $htmlFh "</table>\n";


    # start error by position
    #print $htmlFh "<a name=\"ErrPos\"><h3>Error By Position</h3></a>\n";
    #print $htmlFh "<h4>Summary</h4><br>\n";
    #$plotFile = sprintf("position_error_summary.txt");
    #&writeTextFile($htmlFh, $plotFile) if (-e $plotFile);

    # start off target analysis
    print $htmlFh "<a name=\"OffTarget\"><h3>Off Target Analysis</h3></a>\n";
    my $offTargetBed = $analysisName . ".offtarget.bed";
    my $offTargetBedSize = 0;
    $offTargetBedSize = -s $offTargetBed;
    if (!$offTargetBedSize || $offTargetBedSize == 0) {
	print $htmlFh "Off target analysis produced no enrichment region, with min_depth=100, min_size=10)<br>\n";
    } else {
	my $offTargetBam = $analysisName . ".offtarget.bam";
	my $offTargetBai = $offTargetBam . ".bai";
	print $htmlFh "<a href=\"$offTargetBam\">Off target bam file</a> (off-target-read bam file, min_depth=100, min_size=10)<br>\n";
	print $htmlFh "<a href=\"$offTargetBai\">Off target bam index file</a><br>\n";
	print $htmlFh "<a href=\"$offTargetBed\">Off target bed</a> (off-target-read enriched regions)<br>\n";
	print $htmlFh "<br>";
	my $offTargetProxBed = $analysisName . ".offtarget.proximity.bed";
	print $htmlFh "<a href=\"$offTargetProxBed\">Off target proximity bed and possible mis-annealing region</a> (off-target-read enriched regions with nearby target region)<br>\n";
	print $htmlFh "<center>\n";
	my $offTargetSnapTab = $analysisName . ".offtarget.snap.tab";
	$plotFile = $offTargetSnapTab;
	#$plotFile = sprintf("offtarget.snap.tab");
	print $htmlFh "<a href=\"$plotFile\">Download</a><br>\n" if(-e $plotFile);
	print $htmlFh "</center>\n";
        print $htmlFh "<right> Click column headers to sort </right><br>";
	print $htmlFh "<center>\n";
	#&writeTabFile($htmlFh, $plotFile) if (-e $plotFile);
	&writeTabFile($htmlFh, $plotFile, 1) if (-e $plotFile);
	print $htmlFh "</center>\n";


    }

    # start immediate 5' base analysis
    print $htmlFh "<a name=\"ImmediateBase\"><h3>Immediate 5-prime Base Analysis</h3></a>\n";
    print $htmlFh "The analysis of immediate 5-prime base from reference which is right next to the start position of the aligned read.<br>\n";
    print $htmlFh "<center>\n";
    print $htmlFh "<table border=1 cellpadding=3 width=50% align='center'>\n";
    print $htmlFh "    <tr>\n";
    print $htmlFh "<td>";
    print $htmlFh "<a href=\"$plotFile\">On-target reads</a><br>\n" if(-e $plotFile);
    #$plotFile = $analysisName . ".ontarget.first.5.prime.base.tab";
    $plotFile = $analysisName . ".ontarget.first.5.prime.base.tab";
    &writeTabFile($htmlFh, $plotFile) if (-e $plotFile);
    print $htmlFh "</td>\n" if(-e $plotFile);
    print $htmlFh "<td>";
    #$plotFile = $analysisName . ".offtarget.first.5.prime.base.tab";
    $plotFile = $analysisName . ".offtarget.first.5.prime.base.tab";
    print $htmlFh "<a href=\"$plotFile\">Off-target reads</a><br>\n" if(-e $plotFile);
    &writeTabFile($htmlFh, $plotFile) if (-e $plotFile);
    print $htmlFh "</td>\n" if(-e $plotFile);
    print $htmlFh "    </tr>\n";
    print $htmlFh "</table>\n";
    print $htmlFh "</center>\n";


    finishHtml($htmlFh);
    chdir $cwd;
    return;
}

sub finishHtml {
    my $fh = shift(@_);
    #print $fh "<br/>\n<form id='exportform' method='POST' action='export/export.php'>\n";
    #print $fh "<input id='exportdata' name='exportdata' type='hidden'/>\n";
    #print $fh "<input id='exportfn' name='exportfn' type='hidden' value='whatever.csv'/>\n</form>\n";
    #open (HTMLP2, "$plugindir/export_csv_html_p2.txt") || die "Can't open $plugindir/export_csv_html_p2.txt: $!";
    #while (<HTMLP2>) {
	#print $fh $_;
    #}
    #close (HTMLP2);

    print $fh "</body>\n";
    print $fh "</html>\n";
    close $fh;
}

sub writeHtmlHeader() {
    my $title = shift(@_);
    my $header = shift(@_);
    my $fh = shift(@_);
    print $fh "<html>\n";
    print $fh "<head>\n";
    print $fh "<title>$title</title>\n";

    #table-sort/export part
    open (HTMLP0, "$plugindir/export_csv_html_p0.txt") || die "Can't open $plugindir/export_csv_html_p0.txt: $!";
    while (my $line = <HTMLP0>) {
	print $fh $line;
    }
    close (HTMLP0);

    print $fh "</head>\n";
    print $fh "<body>\n";
    print $fh "<h1><center>$header</center></h1>\n";
}

sub runAndCheckCommand {
    if ($justhtml) {
	$retVal = 0;
	return;
    }

    my $con = shift(@_);
    print($command, "\n");
    $retVal += system($command);
    print "At ", `date`, "return value now is $retVal\n";
    if ($retVal != 0) {
	if ($con) {
	    $retVal = 0;
	    print "Program encountered problem, but keep going.\n";
	} else {
	    print "Program stop after running.\n";
	    exit(0);
	}
    }
}

sub getPriorityData {
    my ($bam, $hotspot, $extbed, $bed, $cut, $p1, $p9) = @_;
    my %pri;
    my %depth;
    my %nonpass;
    my %pass;
    my %coor;

    open (EXT, $extbed) || die "Can't open $extbed: $!\n";
    while (<EXT>) {
        chomp;
        my @d = split(/\t/, $_);
        $pri{$d[0]} = $d[3];
    }
    close (EXT);
    open (BED, $bed) || die "Can't open $bed: $!\n";
    while (<BED>) {
        chomp;
        my @d = split(/\t/, $_);
        $coor{$d[3]}{chr} = $d[0];
        $coor{$d[3]}{start} = $d[1];
        $coor{$d[3]}{stop} = $d[2];
    }
    close (BED);
    my $depth = "$samtools depth -b $hotspot $bam |";
    open (DEP, $depth) || die "Can't open samtools depth output:\n\t-- $depth\n";
    my %depthCt;
    while (<DEP>) {
        chomp;
        my @d = split(/\t/, $_);
        $depth{$d[0]}{$d[1]} = $d[2];
    }
    close (DEP);


    my %dupCheck;
    open (HOT, $hotspot) || die "Can't open $hotspot: $!\n";
    while (<HOT>) {
        next if (/^track/);
        chomp;
        my @d = split(/\t/, $_);
        my $chr = $d[0];
        my $start = $d[1];
        my $stop = $d[2];
	if ($dupCheck{$start . '_' . $stop}) {
	    #print "$chr\t$start\t$stop\tduplicated\n";
	    next;
	}
	$dupCheck{$start . '_' . $stop} = 1;
        if ($start == $stop) { #insertion notations, not listed in extended.bed file anyway
            #print "$chr\t$start\t$stop\tinsertion in bed file\n";
            next;
        }
        my $avgDepth = 0;
        my $totalDepth = 0;
        foreach my $pos (keys %{ $depth{$chr} }) {
            $totalDepth += $depth{$chr}{$pos} if ($pos >= ($start + 1) && $pos <= $stop);
        }
        $avgDepth = $totalDepth / ($stop - $start);

        foreach my $amp (keys %coor) {
            if ($coor{$amp}{chr} eq $chr && ($coor{$amp}{start} <= $start && $coor{$amp}{stop} >= $stop)) {
                if ($avgDepth >= $cut) {
                    $pass{$pri{$amp}}++;
                } else {
                    $nonpass{$pri{$amp}}++;
                }
            }
        }
    }
    close (HOT);

    $pass{1} = $pass{1} ? $pass{1} : 0;
    $pass{9} = $pass{9} ? $pass{9} : 0;
    $nonpass{1} = $nonpass{1} ? $nonpass{1} : 0;
    $nonpass{9} = $nonpass{9} ? $nonpass{9} : 0;
    print "pass/nonpass pri 1 is $pass{1}, $nonpass{1} and pri 9 is $pass{9}, $nonpass{9}\n";

    open (P1, ">$p1") || die "Can't open $p1: $!\n";
    open (P9, ">$p9") || die "Can't open $p9: $!\n";
    foreach my $k (1, 9) {
        my $perc = sprintf("%.3f", (100 * $pass{$k}/ ($pass{$k} + $nonpass{$k})));
        print P1 "$perc\n" if ($k == 1);
        print P9 "$perc\n" if ($k == 9);
    }
    close (P1);
    close (P9);
}

sub findPreviousOutputFolder {
    my $currFolder = shift;
    my $preNum = shift;
    my $upperDir = $currFolder . '/../';
    my $preOutDir = $currFolder;
    $preOutDir =~ s/\_out\.\d+/\_out\.$preNum/;
    return $preOutDir;
}
