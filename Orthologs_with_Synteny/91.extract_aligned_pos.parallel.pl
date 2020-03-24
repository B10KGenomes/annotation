#!/usr/bin/perl -w

=head1 Description

This script is used to create the work shell which is used to extract positions from HAL.

=head1 Contact & Version

  Author: Qi Fang, fangqi@genomics.cn
  Version: 1.0,  Date: 2019-01-17
  Version: 1.1,  Date: 2019-02-18, 1. --unique of hal2maf will lost some columns, deleted; 2. add 03.splitLocallyColinear.pl to split some inversions on the same scaffold
  Version: 2.0,  Date: 2019-09-11 Modified: discard Method 1 and 2; add --ortholog to output the longest and most identical locally colinear alignments and --species_list
  Version: 2.1,  Date: 2019-10-21 Modified: add hal2maf detection form environment
  Version: p3,   Date: 2019-11-19 Modified from extract_shell.pl

=head1 Command-line Option

  --hal <str>       HAL alignments file
  --input <str>     input genome length BED format file
  --ref <str>       reference of multiple alignments
  --tar <str>       target of multiple alignments
  --split_dir <str> output directory

  --help            output help information to screen

=head1 Usage Exmples

  perl extract_shell.pl -i Anc001.fasta.len.bed -r Anc001 --split_dir Anc001_extract --hal allspecies.hal

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin);

my ($Input,$Ref,$Tar,$Split_dir,$Hal,$Help);
GetOptions(
        "ref:s"=>\$Ref,
        "tar:s"=>\$Tar,
        "input:s"=>\$Input,
        "split_dir:s"=>\$Split_dir,
        "hal:s"=>\$Hal,
        "help"=>\$Help
);

die `pod2text $0` if ($Help);

open IN,"<$Input";
mkdir "$Split_dir" unless (-d "$Split_dir");
my $parentdir="00";
my $subdir = "000";
my $parentloop=0;
my $loop = 0;
my $cmd;

while(<IN>){
    if($loop % 100 == 0){
        if($parentloop % 100 ==0){
            $parentdir++;
            mkdir ("$Split_dir/$parentdir");
            $subdir="000";
        }
        $subdir++;
        mkdir("$Split_dir/$parentdir/$subdir");
        $parentloop++;
    }
    chomp;
    my @a = split(/\s+/,$_);
    my $qr_file = "$Split_dir/$parentdir/$subdir/$a[3].bed";
    open OUT, ">$qr_file" || die "fail creat $qr_file";
    print OUT "$a[0]\t$a[1]\t$a[2]\t$a[3]\n";
    close OUT;
    my $HalLiftover = `which halLiftover`;
    chomp $HalLiftover;
    $cmd .= "$HalLiftover $Hal $Ref $qr_file $Tar $qr_file.$Tar.bed\n";
    $loop++;
}

open OUT, ">split_halLiftover_$Ref.sh" || die "fail creat split_halLiftover_$Ref.sh";
print OUT $cmd;
close OUT;
