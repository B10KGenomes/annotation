#!/usr/bin/perl -w

=head1 Description

This script is used to identify the gene order rank, and deal with the duplicates.

=head1 Contact & Version

  Author: Yuan Deng & Qi Fang, dengyuan@genomics.cn & fangqi@genomics.cn
  Version: p1,   Date: 2019-11-19
  Version: p2,   Date: 2020-01-12 Modified by Qi Fang, to differentiate the synteny support into perfect and good.

=head1 Usage Exmples

  perl 95.syntenic_ort.pl intersect.ortholog.list.order reference.gff query.gff > intersect.ortholog.list.order.ort

=cut

use strict;
use File::Basename;

die "Usage: <IN> <ref.gff> <target.gff>\n" unless @ARGV == 3;

## read gff
my %chrGeneNum;
my @gff_files=($ARGV[1],$ARGV[2]);
foreach my$file(@gff_files){
    open IN, $file;
    while (<IN>) {
        my @info = split /\s+/;
        next unless $info[2] eq "mRNA";
        my $sp = $1, if($info[8]=~/ID=([A-Z]+)_\w+;/);
        $chrGeneNum{$sp}{$info[0]} ++;
    }
    close IN;
}

## restore the whole table
my @data;
open IN, $ARGV[0];
@data = <IN>;
close IN;
## creat psudo lines
my @tmp = split /\s+/, $data[0];
my $col_n = @tmp;
@tmp = ();
for (my $i = 0; $i < $col_n; $i ++) {
    push @tmp, 10e10;
}
my $tmp = join "\t", @tmp;
## add five pseudo lines to the head and the tail
for (my $i = 0; $i < 5; $i ++) {
    unshift @data, "$tmp\n";
    push @data, "$tmp\n";
}
for (my $i = 5; $i < @data - 5;  $i ++) {
    my @info = split /\s+/, $data[$i];
    next if $info[7] eq "NA";

    my $up_check = 0;
    my $up_i;
    for (my $j = $i - 1; $j >= 0; $j --) {
        my @up = split /\s+/, $data[$j];
        next if $up[7] eq "NA";
        next, if($info[0] eq $up[0]); ##add by fangqi, deal with #many# to many in cactus, same ref gene
        next, if($info[6] eq $up[6]); ##add by fangqi, deal with #many# to many in cactus, same ref chr to same gene
        if($up[7] ne $info[7]){ ##add by fangqi, deal with many to #many# in cactus, same ref chr to diff query chr
            next, if($j-1<0);
            my @up_tmp = split /\s+/, $data[$j-1];
            if($up_tmp[0] eq $up[0]){next;}
        }
        if ($up[1] eq $info[1] && $up[7] eq $info[7] ) {
            if ( abs($up[8] - $info[8]) == 1) {  ##add by fangqi, deal with nearly insertion
                $up_check = 2; ##add by fangqi, deal with nearly insertion
            } elsif ( abs($up[8] - $info[8]) < 5) {
                $up_check = 1;
            }
            $up_i = $j;
            last;
        } else {
            last;
        }
    }

    my $down_check = 0;
    my $down_i;
    for (my $j = $i + 1; $j < @data; $j ++) {
        my @down = split /\s+/, $data[$j];
        next if $down[7] eq "NA";
        next, if($info[0] eq $down[0]); ##add by fangqi, deal with #many# to many in cactus, same ref gene
        next, if($info[6] eq $down[6]); ##add by fangqi, deal with #many# to many in cactus, same ref chr to same gene
        if($down[7] ne $info[7]){ ##add by fangqi, deal with many to #many# in cactus, same ref chr to diff query chr
            next, if($j+1>$#data);
            my @down_tmp = split /\s+/, $data[$j+1];
            if($down_tmp[0] eq $down[0]){next;}
        }
        if ($down[1] eq $info[1] && $down[7] eq $info[7] ) { 
            if ( abs($down[8] - $info[8]) == 1 ) {  ##add by fangqi, deal with nearly insertion
                $down_check = 2;  ##add by fangqi, deal with nearly insertion
            } elsif ( abs($down[8] - $info[8]) < 5 ) {
                $down_check = 1;
            }
            $down_i = $j;
            last;
        } else {
            last;
        }
    }
    
    my $chr_gene_num_check = 0;
    for (my $j = 0; $j < @info-1; $j += 6) {
        my $sp = (split /_/, $info[$j])[0];
        my $chr = $info[$j+1];
        $chr_gene_num_check = 1 if $chrGeneNum{$sp}{$chr} == 1;
    }
    
    if ($up_check && $down_check) {
        chomp($data[$i]);
        if ($up_check==2 && $down_check==2){
            print "$data[$i]Perfect\t\n";
        } elsif ($up_check==2 or $down_check==2) {
            print "$data[$i]Good\t\n";
        } else {
        #print $data[$up_i], "#$data[$i]", $data[$down_i], "\n";
            print "$data[$i]OK\t\n";
        }
    } elsif ($up_check) {
        #print $data[$up_i], "#$data[$i]", "\n";
        chomp($data[$i]);
        print "$data[$i]DownNone\t\n";
    } elsif ($down_check) {
        #print "#$data[$i]", $data[$down_i], "\n";
        chomp($data[$i]);
        print "$data[$i]UpNone\t\n";
    } elsif ($chr_gene_num_check) {
        chomp($data[$i]);
#        my @tmp = split /\s+/, $data[$i];
#        if($tmp[8]==1){
            print "$data[$i]Orphan\t\n";
#        }else{
#            print "$data[$i]Indel\t\n";
#        }
    }
}
