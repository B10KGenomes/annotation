#!/usr/bin/perl -w
use strict;

if (@ARGV != 2)
{
	print "perl $0 <.gff> <ident cutoff>\n";
	exit;
}

use strict;

my $gff = shift;
my $cut= shift;
my %hshp;

open IN,$gff or die;
while (<IN>)
{
	my @a = split;
	if ($a[2] eq "mRNA"){
		next if ($a[5] < $cut);
		/ID=([^;]+);/;
		$hshp{$1} = $a[5];
		print $_;
	}elsif($a[2] eq "CDS"){
		/Parent=([^;]+);/;
		next unless ($hshp{$1});
		print $_;
	}
}
close IN;
