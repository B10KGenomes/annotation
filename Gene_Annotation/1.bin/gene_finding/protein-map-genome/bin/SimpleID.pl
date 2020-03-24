#!/usr/bin/perl 
use strict;


my $infile=shift;

my $output;
open IN,$infile or die "Fail $infile:$!";
$/=">";<IN>;$/="\n";
while(<IN>){
	my ($id,$seq);
	if (/^(\S+)/){
		$id=$1;
		if ($id=~/\|/){
			die "Found \"|\" in ID:line $. .\n\n";
		}
	}else{
		die "Format ERROR in line $. .\n\n";
	}
	$/=">";
	$seq=<IN>;
	chomp $seq;
	$output.=">$id\n$seq";
	$/="\n"; 
} 
close IN;

my $outfile=$infile.".SimpleID.fa";
open OUT,">$outfile" or die "Fail $outfile:$!";
print OUT $output;
close OUT;
