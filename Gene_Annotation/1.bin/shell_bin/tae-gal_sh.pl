#!/usr/bin/perl

#2011.10.03

if (@ARGV != 2)
{
	print "perl $0 <dir> <gal or tae>\n";
	exit;
}

use strict;
use FindBin qw($Bin $Script);

my $dir = shift;
my $mark = shift;

`rm  $dir/*genewise.gff.iden*`;
`rm $dir/*filter.$mark`;
`rm -r $dir/cluster*`;

my @identitY_pw = <$dir/*.list.filter>;
my @gff = <$dir/*genewise.gff>;
my @cds = <$dir/*cds>;
`perl $Bin/../common_bin/pseudogene_filter.pl $gff[0] $cds[0] > $gff[0].PF.gff 2> $gff[0].pseudo.gff`;

`perl $Bin/../common_bin/algnRat2identity_gff.pl $identitY_pw[0] $gff[0].PF.gff >  $gff[0].PF.ident.gff`;

`perl $Bin/../common_bin/fishInWinter.pl $Bin/../../0.homolog_data/gal_tae.ref.orth.$mark $identitY_pw[0] -fc 2 > $identitY_pw[0].orth`;
`perl $Bin/../common_bin/fishInWinter.pl $Bin/../../0.homolog_data/gal_tae.ref.para.$mark $identitY_pw[0] -fc 2 > $identitY_pw[0].para`;

`perl $Bin/../common_bin/fishInWinter.pl -ff gff $identitY_pw[0].orth $gff[0].PF.ident.gff > $gff[0].PF.ident.orth.gff`;
`perl $Bin/../common_bin/fishInWinter.pl -ff gff $identitY_pw[0].para $gff[0].PF.ident.gff > $gff[0].PF.ident.para.gff`;


