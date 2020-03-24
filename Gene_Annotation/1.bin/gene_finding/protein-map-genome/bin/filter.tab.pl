#!/usr/bin/perl

use strict;

die "perl $0 <tab> [value] > outfile \n" if @ARGV < 1;

my $Tab=shift;
my $cutoff=shift;

$cutoff ||= 0.8;

my %tab;
open IN,$Tab or die $!;
while(<IN>){
	next if /^#/;
	chomp;
	my @a=split /\t+/;
	my $id=$a[0];
	$id=$1 if($id =~ /(\S+)-D\d+$/);
	push @{$tab{$id}},[@a];
}
close IN;

for my $id(keys %tab){
	@{$tab{$id}}=sort {$b->[10] <=> $a->[10]} @{$tab{$id}};
	my $p=$tab{$id};
	my $best=$p->[0]->[10];
	for (@{$tab{$id}}){
		print join("\t",@$_)."\n" if($best*$cutoff <= $_->[10]);
	}
}
