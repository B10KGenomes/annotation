#!/usr/bin/perl -w
use strict;
#use Getopt::Long;
#use File::Basename;
sub usage{
    print STDERR <<USAGE;
    ##############################################
            Version 1.0 by Wing-L  2011.12.28

      usage: $0 <gff> (<print>)>STDERR(average stat) STDOUT(every genes)

        default,print out the average statistics of gff in STDERR,each columns are
                1        2        3            4               5             6               7                 8
            Gene_num  GeneLen  ExonNum  GeneExonTotalLen  SingleExonLen  IntronNum  GeneIntronTotalLen  SingleIntronLen

      if set <print>,will print information of each genes in STDOUT.And the first colum is
          GeneID,othes as default.
    ##############################################
USAGE
    exit;
}
&usage if(@ARGV <1 or @ARGV>2);
my ($gff,$print)=@ARGV;

############# global variable
#my ();

#GetOptions(
#  ""=>\,
#  ""=>\,
#  ""=>\,
#  ""=>\,
#);
############################  BEGIN MAIN  ############################
my %stat;
open IN,"<$gff" or die("$!:$gff\n");
while(my $line=<IN>){
    my @info=(split /\s+/,$line);
    next if($line=~/\#/ or !defined $info[1]);
    $info[8]=~/=([^;]+);/;
    if($info[2] eq 'mRNA'){
        $stat{$1}{'mRNA'}+=$info[4]-$info[3]+1;
    }elsif($info[2] eq 'CDS'){
        ++$stat{$1}{'CDS'}[0];
        $stat{$1}{'CDS'}[1]+=$info[4]-$info[3]+1;
    }
}
close IN;
my $num=keys %stat;
my ($total_exon_num,$total_exon_len,$total_intron_num,$total_intron_len,$total_len);
foreach my $e (keys %stat){
    my $intron_num=$stat{$e}{'CDS'}[0]-1;
    my $intron_len=$stat{$e}{'mRNA'}-$stat{$e}{'CDS'}[1];
    my $average_intron_len=$intron_num>0?$intron_len/$intron_num:0;
    my $average_exon_len=$stat{$e}{'CDS'}[1]/$stat{$e}{'CDS'}[0];
    printf "$e\t$stat{$e}{'mRNA'}\t$stat{$e}{'CDS'}[0]\t$stat{$e}{'CDS'}[1]\t%.2f\t$intron_num\t$intron_len\t%.2f\n",$average_exon_len,$average_intron_len if(defined $print);
    $total_len+=$stat{$e}{'mRNA'};
    $total_exon_num+=$stat{$e}{'CDS'}[0];
    $total_exon_len+=$stat{$e}{'CDS'}[1];
    $total_intron_num+=$intron_num;
    $total_intron_len+=$intron_len;
}
my $average_len=$total_len/$num;
my $average_exon_num=$total_exon_num/$num;
my $gene_average_exon_len=$total_exon_len/$num;
my $single_average_exon_len=$total_exon_len/$total_exon_num;
my $average_intron_num=$total_intron_num/$num;
my $gene_average_intron_len=$total_intron_len/$num;
my $single_average_intron_num=$total_intron_len/$total_intron_num;
printf STDERR "$num\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",$average_len,$average_exon_num,$gene_average_exon_len,$single_average_exon_len,$average_intron_num,$gene_average_intron_len,$single_average_intron_num;
############################   END  MAIN  ############################

###############  sub  ###############
