#!/usr/bin/perl -w

=head1 Description

This script is used to get all gene pairs and thier overlaped length.

=head1 Contact & Version

  Author: Qi Fang, fangqi@genomics.cn
  Version: p1,   Date: 2019-11-19

=head1 Usage Exmples

  perl 92.get_overlap.pl pairs.list reference.cds.fa query.cds.fa > results

=cut

use strict;

my $input = shift;
my $ref_cds = shift;
my $tar_cds = shift;

open IN,$tar_cds;
my %cds_len;
while(<IN>){
    chomp;
    my @a = split /\t+/;
    my $id = $a[3];
    $cds_len{$id}+=$a[2]-$a[1];
}close IN;

open IN,$ref_cds;
while(<IN>){
    chomp;
    my @a = split /\t+/;
    my $id = $a[3];
    $cds_len{$id}+=$a[2]-$a[1];
}close IN;

open IN,$input;
my %int_len;
while(<IN>){
    chomp;
    my @a = split /\t+/;
    my $id2 = $a[3];
    my $id1 = $a[7];
    $int_len{$id1}{$id2}+=$a[-1];
}close IN;

foreach my $id1 (sort keys %int_len){
    foreach my $id2 (sort keys %{$int_len{$id1}}){
        my $re = $int_len{$id1}{$id2} / (($cds_len{$id1}+$cds_len{$id2})/2);
#        next, if($int_len{$id1}{$id2}<0.03);
        print "$id1\t$id2\t$re\n";
    }
}

####################################################
################### Sub Routines ###################
####################################################

##read sequences in fasta format and calculate length of these sequences.
sub read_fasta{
    my ($file,$p)=@_;
    if($file =~ /\.gz$/){
        open IN,"gzip -dc $file|" or die $!;
    }else{
        open IN,$file or die "Fail $file:$!";
    }
    $/=">";<IN>;$/="\n";
    while(<IN>){
        my ($id,$seq);
        if (/^(\S+)/){
            $id=$1;
        }
        $/=">";
        $seq=<IN>;
        chomp $seq;
        $seq=~s/\s//g;
        $p->{$id}=length($seq);
        $/="\n";
    }
    close IN;
}

sub read_gff{
    my ($file,$p,$chr,$strand)=@_;
    open IN,$file or die "Fail $file:$!";
    while(<IN>){
        chomp;
        my @a = split /\s+/;
        next, unless($a[2] eq "CDS");
        my $name = $1, if($a[8]=~/Parent=([A-Z]+_R[0-9]+);Source=.*/);
        $strand->{$name}=$a[6];
        push @{$p->{$name}}, ($a[3]-1,$a[4]);
        push @{$chr->{$a[0]}}, $name;
    }
    close IN;
}

sub read_rna_len{
    my ($file,$rna_len)=@_;
    open IN,$file or die "Fail $file:$!";
    while(<IN>){
        chomp;
        my @a = split /\s+/;
        my $name = $1, if($a[3]=~/ID=([A-Z]+_R[0-9]+);Source=.*/);
        $rna_len->{$name}=($a[2]-$a[1]);
    }
    close IN;
}
