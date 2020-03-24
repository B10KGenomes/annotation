#!/usr/bin/perl
if (@ARGV != 1)
{
        print "perl $0 <human_dir>\n";
        exit;
}
use strict;
use FindBin qw($Bin $Script);
my $dir=shift;
`cd $dir`;

### filter unqualified UniForm_genes
`cat $dir/*.pep.*.gff.ident40.PF.gff |  awk '$3=="mRNA"{print $NF}' | sed 's/ID=//;s/;.*//' | awk -F"-" 'BEGIN{while (getline < "$Bin/uniform_gene_set_exTRAN.list"){da[$1]=$1}}{print $0"\t"da[$1]}' | awk '$2!=""{print $1}' | perl $Bin/../common_bin/fishInWinter.pl -bf table -ff gff - $dir/*.pep.*.gff.ident40.PF.gff > human.ident40.PF.UF.gff`;

### filter Multiple2Single
`cat human.ident40.PF.UF.gff | awk '$3=="CDS"{print $NF}' | sed 's/Parent=//;s/;$//' | uniq -c | awk '$1==1{print $2}' | awk -F"-" 'BEGIN{while (getline < "$Bin/uniform_gene_set_exTRAN.or.multipleExon.list"){da[$1]=$1}}{print $0"\t"da[$1]}' | awk '$2!=""{print $1}' | perl $Bin/../common_bin/fishInWinter.pl -bf table -ff gff -except -  human.ident40.PF.UF.gff >  human.ident40.PF.UF.M2S.gff`;

### filter mRNA shorter than 150
`cat  human.ident40.PF.UF.M2S.gff | awk '$3=="mRNA"{x=$5-$4; print $NF"\t"x}' | sed 's/ID=//;s/;.*\t/\t/' | awk '$2>=150{print $1}' | perl $Bin/../common_bin/fishInWinter.pl -bf table -ff gff -  human.ident40.PF.UF.M2S.gff >  human.ident40.PF.UF.M2S.len150.gff`;

### hcluster for tae_gal.gff
`perl $Bin/clustergff/clustergff_score_b10k.pl $dir/human.ident40.PF.UF.M2S.len150.gff --run qsub --cutoff 0.4 --cpu 50`;

`cd -`;
