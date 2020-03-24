#!/usr/bin/perl

if (@ARGV != 5)
{
        print "perl $0 <homo_gff> <human_gff> <trans_gff> <repeatmasked_genome> <Gene_Id>\n";
        exit;
}
use strict;
use FindBin qw($Bin $Script);

my $gff1=shift; ##homo
my $gff2=shift; ##human
my $gff3=shift; ##trans
my $fasta=shift; ##masked
my $id=shift;

`ln -s $gff1 homo.gff`;
`ln -s $gff2 human.gff`;
`ln -s $gff3 trans.gff`;

### get gene list
`cat homo.gff | awk '$3=="mRNA"' | sed 's/.*ID=//;s/;.*//' > homo.gff.list`;
`cat human.gff | awk '$3=="mRNA"' | sed 's/.*ID=//;s/;.*//' > human.gff.list`;
`cat trans.gff | awk '$3=="mRNA"' | sed 's/.*ID=//;s/;.*//' > trans.gff.list`;
### dup10
`cat homo.gff.list | sed 's/-D.*//' | sort | uniq -c | awk '$1>10{print $2}' | grep -f - homo.gff.list > homo.gff.list.dup10`;
`cat human.gff.list | sed 's/-D.*//' | sort | uniq -c | awk '$1>10{print $2}' | grep -f - human.gff.list > human.gff.list.dup10`;
`cat trans.gff.list | sed 's/-D.*//' | sort | uniq -c | awk '$1>10{print $2}' | grep -f - trans.gff.list > trans.gff.list.dup10`;
### singleExon
`perl $Bin/get_single_exon_gene.pl homo.gff > homo.gff.list.singleExon`;
`perl $Bin/get_single_exon_gene.pl human.gff > human.gff.list.singleExon`;
`perl $Bin/get_single_exon_gene.pl trans.gff > trans.gff.list.singleExon`;
### repeat
`perl $Bin/../common_bin/getGene.pl homo.gff $fasta > homo.gff.cds`;
`perl $Bin/../common_bin/getGene.pl human.gff $fasta > human.gff.cds`;
`perl $Bin/../common_bin/getGene.pl trans.gff $fasta > trans.gff.cds`;
`perl $Bin/masked_rate.pl homo.gff.cds | awk '$2>0.7{print $1}' > homo.gff.cds.repeat70;`;
`perl $Bin/masked_rate.pl human.gff.cds | awk '$2>0.7{print $1}' > human.gff.cds.repeat70;`;
`perl $Bin/masked_rate.pl trans.gff.cds | awk '$2>0.7{print $1}' > trans.gff.cds.repeat70;`;
### filter
`perl $Bin/get_filter_list.pl homo.gff.list.dup10 homo.gff.list.singleExon homo.gff.cds.repeat70 > homo.gff.list.removed`;
`perl $Bin/get_filter_list.pl human.gff.list.dup10 human.gff.list.singleExon human.gff.cds.repeat70 > human.gff.list.removed`;
`perl $Bin/get_filter_list.pl trans.gff.list.dup10 trans.gff.list.singleExon trans.gff.cds.repeat70 > trans.gff.list.removed`;
### rest gff
`perl $Bin/../common_bin/fishInWinter.pl -bf table -ff gff -except homo.gff.list.removed homo.gff > homo.gff.rest.gff`;
`perl $Bin/../common_bin/fishInWinter.pl -bf table -ff gff -except human.gff.list.removed human.gff > human.gff.rest.gff`;
`perl $Bin/../common_bin/fishInWinter.pl -bf table -ff gff -except trans.gff.list.removed trans.gff > trans.gff.rest.gff`;
### add without overlaps
`perl $Bin/FindOverlapAtCDSlevel.pl homo.gff.rest.gff human.gff.rest.gff | sed '1d' | awk '{print $2}' | perl $Bin/../common_bin/fishInWinter.pl -bf table -ff gff -except - human.gff.rest.gff | cat homo.gff.rest.gff - > ${id}.tmp1.gff`;
`perl $Bin/FindOverlapAtCDSlevel.pl ${id}.tmp1.gff trans.gff.rest.gff | sed '1d' | awk '{print $2}' | perl $Bin/../common_bin/fishInWinter.pl -bf table -ff gff -except - trans.gff.rest.gff | cat ${id}.tmp1.gff - | perl $Bin/get_taxon.pl - > ${id}.homolog.gff`;
