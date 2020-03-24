#!usr/bin/perl -w
use strict;

my $dup10=shift;
my $sExon=shift;
my $repeat=shift;
my %evidence;

open (IN1,"$sExon") or die $!;
while (<IN1>){
	chomp;
	$evidence{$_}=$_;
}
close IN1;

open (IN2,"$repeat") or die $!;
while (<IN2>){
	chomp;
	$evidence{$_}=$_;
}
close IN2;

open (IN3,"$dup10") or die $!;
while (<IN3>){
	chomp;
	print "$_\n" if (exists $evidence{$_});
}
close IN3;
