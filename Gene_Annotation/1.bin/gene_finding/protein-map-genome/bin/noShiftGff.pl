#!/usr/bin/perl

use strict;

my $inlist=shift;
my $gff=shift;

my %list;
open IN,$inlist or die $!;
while(<IN>){
	next if /^#/;
	chomp;
	$list{$_}=1;
}
close IN;

my $num=0;
open IN,$gff or die $!;
while(<IN>){
	chomp;
	my $line=$_;
	my @a=split /\t/,$line;
	my $id=$1 if($a[8] =~ /ID=([^;\s]+)/ && $a[2] eq 'mRNA');
	if($a[2] eq 'mRNA' && exists $list{$id}){
		$a[8] =~ s/(\S+Shift=)(\d+)(\S+)/$1$num$3/;
	}
	print join("\t",@a)."\n";
}
close IN;
