#!/usr/bin/perl -w
use strict;

open (IN,"$ARGV[0]") or die $!;
while (<IN>){
	my @t = split;
	if ($t[-1] =~ /ENSTGUP/ || $t[-1] =~ /ENSGALP/){
		$t[1] = "GeneWise";
	}elsif($t[-1] =~ /ENSP00/){
		$t[1] = "HumanGene";
	}else{
		$t[1] = "TRANS";
	}
	my $l = join"\t",@t;
	print "$l\n";
}
close IN;
