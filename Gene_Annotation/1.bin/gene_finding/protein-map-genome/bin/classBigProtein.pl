#!/usr/bin/perl

use strict;

die "perl $0 <shell> <pep.fa.len> [outdir]\n" if @ARGV < 2;

my $Shell=shift;
my $Len=shift;
my $Outdir=shift;


$Outdir ||= ".";
$Outdir =~ s/\/$//;

my %len;
open IN,$Len or die $!;
while(<IN>){
	next if /^#/;
	$len{$1}=$2 if(/(\S+)\s+(\S+)/);
}
close IN;

my $name = $1 if( $Shell =~ /([^\/]+)$/);
open IN,$Shell or die $!;
open OUT1,">$Outdir/$name.bt1k.shell" or die $!;
open OUT2,">$Outdir/$name.st1k.shell" or die $!;
open OUT3,">$Outdir/$name.bt3k.shell" or die $!;
while(<IN>){
	next if /^#/;
	my @a=split /\s+/;
	my $id;
	$id = $1 if ($a[7] =~ /([^\/]+)\.fa$/); ##$a[5] --> $a[7],because of the genewise libpath
	$id = $1 if ($a[7] =~ /([^\/]+)-D\d+\.fa$/); ##$a[5] --> $a[7]
	if($len{$id} >= 3000){
		print OUT3 $_;
	}elsif($len{$id} >= 1000){
		print OUT1 $_;
	}else{
		print OUT2 $_;
	}
}
close IN;
close OUT1;
close OUT2;
close OUT3;
