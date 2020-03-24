#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(dirname basename);
use FindBin qw($Bin $Script);

=head1 Usage
 perl clustergff_score_b10k.pl [options] <gff>
 --cpu <int>          set the cpu number to use in parallel, default=8
 --run <str>          set the parallel type, qsub, or multi, default=multi
 --cutoff <int>       set CDSoverlap rate cutoff, default = 0.4
 --verbose            output verbose information to screen,default no
 --help               output help information to screen,default no

=cut

my ($cpu,$run,$cutoff,$verbose,$help);
GetOptions(
        "cpu:i"=>\$cpu,
		"run:s"=>\$run,
		"cutoff:f"=>\$cutoff,
		"verbose!"=>\$verbose,
		"help!"=>\$help
);

$cpu ||= 8;
$cutoff ||= 0.4;
$run ||= "multi";

die `pod2text $0` if (@ARGV < 1 || $help);

my $file = shift;
my $file2 = "$file.temp.gff";

open (IN, $file) || die $!;
open (OUT, ">$file2") || die $!;
while (<IN>){
	chomp;
	next if (/^#/);
	my @c = split /\t/, $_;
	$c[0] = "$c[0]$c[6]";
	my $out = join("\t", @c) . "\n";
	print OUT $out;
}
close IN;
close OUT;

my $basename = basename ($file2, '.gff');
my $dirname = dirname ($file2); 

`perl $Bin/scores_gff.pl $file2 > $dirname/$basename.score`;
`perl $Bin/run_cluster_b10k.pl  $file2 $dirname/$basename.score $run $cutoff $cpu`;

my $file5="$file.uncluster.gff";
open (IN,"$file2.uncluster.gff");
open OUT,">$file5";
while(<IN>){
        chomp;
        my @c = split /\t/, $_;
        $c[0] =~ s/\S$//;
        my $out = join("\t", @c) . "\n";
        print OUT $out;
}
close IN;
close OUT;
`rm $file2.uncluster.gff`;

`mv $file2.dist.hcluster.reference $file.dist.hcluster.reference`;

my $file3 = "$file2.nr.gff";
my $file4 = "$file.nr.gff";
open (IN, $file3) || die $!;
open (OUT, ">$file4") || die $!;
while (<IN>){
	chomp;
	my @c = split /\t/, $_;
	$c[0] =~ s/\S$//;
	my $out = join("\t", @c) . "\n";
	print OUT $out;
}
close IN;
close OUT;

`rm -r $dirname/cluster*` unless ($verbose);
`rm -rf $file.temp.gff*` unless ($verbose);
