#!/usr/bin/perl -w
use strict;
my $input = shift;
my $gff = shift;
my %hash;
open IN,$input;
while(<IN>){
	chomp;
	my @a = split /\t+/;
	$hash{$a[0]}=1;
}close IN;

open IN,$gff;
while(<IN>){
	chomp;
	my @a = split /\t+/;
	my $id;
	if($a[8]=~/ID=(\S+);Shift\S+/){$id = $1;}
	elsif($a[8]=~/Parent=(\S+);/){$id = $1;}
	if(exists $hash{$id}){print "$_\n";}
}close IN;
