#!/usr/bin/perl

if (@ARGV != 1)
{
	print "perl $0 <.gff>\n";
	exit;
}

use strict;

my $fileGff = shift;

&get_scores_gff($fileGff); # sub 1



## sub 1
sub get_scores_gff
{
	my $file = shift;

	open IN,$file or die;
	while (<IN>)
	{
		my @c = split;
		if ($c[2] eq "mRNA")
		{
			$c[8] =~ /ID=([^;]+);/;
			print "$1\t$c[5]\n";
		}
	}
	close IN;
}
