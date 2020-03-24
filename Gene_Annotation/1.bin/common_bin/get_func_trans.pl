#!/usr/bin/perl -w
use strict;

my $func_file = shift;
my $gene_file = shift;

open IN,$func_file;
my %hash;
while(<IN>){
	chomp;
	my @a = split /\s+/;
	$hash{$a[0]} = 1;
}close IN;

open IN,$gene_file;
while(<IN>){
	chomp;
	my @a = split /-/;
	if(exists $hash{$a[0]}){print "$_\n";}
}close IN;
