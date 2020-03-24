#!/usr/bin/perl
use strict;
use warnings;


die "Usage:$0 <net>\n\n" if @ARGV<1;

my $net=shift;

print "#",join("\t",'chr','chr_len','start','end','align_len','strand','chr2','start2','end2','align_len2')."\n";

open IN,$net or die "$!";
my ($chr,$len);
while(<IN>){
	if (/^net\s+(\S+)\s+(\d+)/){
		$chr=$1;
		$len=$2;
	}elsif(/^\sfill\s+(\d+)\s+(\d+)\s+(\S+)\s+([\+-])\s+(\d+)\s+(\d+)/){
		my ($start,$l,$chr2,$str,$start2,$l2)=($1,$2,$3,$4,$5,$6);
		my $end=$start+$l-1;
		my $end2=$start2+$l2-1;
		print join("\t",$chr,$len,$start,$end,$l,$str,$chr2,$start2,$end2,$l2)."\n";
	}
}
close IN;
