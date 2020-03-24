#!/usr/bin/perl
if (@ARGV != 1)
{
        print "perl $0 <trans_dir>\n";
        exit;
}
use strict;
use FindBin qw($Bin $Script);

my $dir=shift;
`cd $dir`;

### remove genes without any functional domain
`perl $Bin/../common_bin/gff_average_len_statistics.pl $dir/TRAN2.pep.*.gff.ident40.PF.gff print | awk '$2>=150 {print $1}' | perl $Bin/../common_bin/get_func_trans.pl $Bin/transcriptome.uniq2.pep.iprscan.gene.ipr_domain_family_only_NoRetro.list - > $dir/trans.ident40.PF.gff.funcRemoved`;
`perl $Bin/../common_bin/get_gff.pl $dir/trans.ident40.PF.gff.funcRemoved $dir/TRAN2.pep.*.gff.ident40.PF.gff > $dir/trans.ident40.PF.len150.funcRemoved.gff`;

### hcluster
`perl $Bin/clustergff/clustergff_score_b10k.pl $dir/trans.ident40.PF.len150.funcRemoved.gff --run qsub --cutoff 0.4 --cpu 50`;

`cd -`;
