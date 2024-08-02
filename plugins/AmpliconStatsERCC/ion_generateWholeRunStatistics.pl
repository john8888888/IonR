#!/usr/bin/perl -w


#read each individual barcode summary_table.txt
#get weighted average Percent no strand bias of all bases
#get weighted Percent end to end read of on target reads


use JSON;


my $currDir = $ARGV[0] || die "Usage: plugin_output_dir\n";

print "add whole run statistics\n";
my @metrics = qw(percent_end_to_end_read_of_on_target_reads percent_no_strand_bias_of_all_bases);

=hold
my $barcodeFile = $currDir . '/../../barcodeList.txt';
my %bc;
open (BF, $barcodeFile) || die "Can't open $barcodeFile: $!\n";
while (<BF>) {
    next if !(/,/);
    my @w = split(/,/, $_);
    $bc{$w[1]}++;
}
close (BF);
=cut

my $jFile = $currDir . '/results.json';
my $jText = from_json(`cat $jFile`);

my %whole;
my %count;
my %digits;

foreach my $k (keys %{ $$jText{'barcodes'} }) {
    #next if !defined($bc{$k});
    foreach my $m (@metrics) {
	next if !defined($$jText{'barcodes'}{$k}{$m});
	my $v = $$jText{'barcodes'}{$k}{$m};
	$v =~ s/\%//;

	my @w = split(/\./, $v);
	$digits{$m} = length($w[1]);

	$v *= 100;
	$v += 0.5;
	$v = int($v);
	$whole{$m} += $v;
	$count{$m}++;
    }
    foreach my $m (@metrics) {
	my $ct = $digits{$m} ? $digits{$m} : 2;
	my $v = $whole{$m} / 100 / $count{$m};
	$$jText{$m} = sprintf "%.${ct}f%%", $v;
    }
}

my $outFile = $currDir . '/results.whole.json';
open (OF, ">$outFile") || die "Can't open $outFile: $!\n";
print OF to_json($jText,{canonical=>1, ascii=>1, pretty=>1});
close (OF);

system("mv $outFile $jFile");
