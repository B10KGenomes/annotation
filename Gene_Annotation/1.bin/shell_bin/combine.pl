#!/usr/bin/perl

if (@ARGV != 3)
{
        print "perl $0 <species_dir> <softmasked.fa> <species_id>\n";
        exit;
}
use strict;
use FindBin qw($Bin $Script);

my $dir=shift;
my $fa=shift;
my $id=shift;

`mkdir $dir/NewFilter`;
`perl $Bin/../NewFilter/NewFilter1.PrimaryFilter.pl $dir/NewFilter`;

`mkdir $dir/combine`;
`cd $dir/combine`;
`perl $Bin/../NewFilter/NewFilter6.Combine.pl $dir/NewFilter/gal_tae.ident40.PF.UF.M2S.len150.gff.nr.gff $dir/HUMAN/human.ident40.PF.UF.M2S.len150.gff.nr.gff $dir/TRAN2/trans.ident40.PF.len150.funcRemoved.gff.nr.gff $fa $id`;
`cd -`;
