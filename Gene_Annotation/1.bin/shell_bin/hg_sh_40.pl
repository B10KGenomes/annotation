#!/usr/bin/perl

#2011.10.03

if (@ARGV != 1)
{
	print "perl $0 <dir>\n";
	exit;
}

use strict;
use FindBin qw($Bin $Script);

my $dir = shift;

my @identitY_pw = <$dir/*.list.filter>;

`awk '\$8 >= 40' $identitY_pw[0] > $identitY_pw[0].40`;

my @gff = <$dir/*genewise.gff>;
`perl $Bin/../common_bin/algnRat2identity_gff.pl $identitY_pw[0].40 $gff[0] >  $gff[0].ident40`;

my @cds = <$dir/*cds>;
`perl $Bin/../common_bin/pseudogene_filter.pl $gff[0].ident40 $cds[0] > $gff[0].ident40.PF.gff 2> $gff[0].ident40.pseudo.gff`;


