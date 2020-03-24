#!/usr/bin/perl -w
use strict;
use File::Basename qw(basename dirname);

die "Usage:perl $0 GFF.file\n" if(@ARGV!=1);

my $GFF_file=shift;
my $name = basename($GFF_file);
my %Gene;
open (IN,$GFF_file) || die "$!\n";
while(<IN>){
	chomp;
	if(/^#/){next;}
	my @c=split /\t/;
	@c[3,4]=@c[4,3] if($c[3]>$c[4]);
	if ($c[2] eq "mRNA"){
		my $id=$1 if($c[8]=~/ID=(\S+?);/ || $c[8]=~/GenePrediction\s+(\S+)/);
		@{$Gene{$id}{mRNA}}=@c;
	}
	if($c[2] eq "CDS"){
		my $id=$1 if($c[8]=~/Parent=(\S+?);/ || $c[8]=~/GenePrediction\s+(\S+)/);
		push @{$Gene{$1}{CDS}},[@c];
	}
}
close IN;

foreach my $id (keys %Gene){
	if (@{$Gene{$id}{CDS}} == 1){print "$id\n"};
}
