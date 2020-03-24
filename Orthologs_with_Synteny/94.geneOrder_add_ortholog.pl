#!/usr/bin/perl

=head1 Description

This script is used to build gene order among all the genes.

=head1 Contact & Version

  Author: Yuan Deng & Qi Fang, dengyuan@genomics.cn & fangqi@genomics.cn
  Version: p1,   Date: 2019-11-19

=head1 Usage Exmples

  perl 94.geneOrder_add_ortholog.pl reference.gff query.gff intersect.ortholog.list > intersect.ortholog.list.order

=cut


if (@ARGV != 3)
{
	print "perl $0 <.gff1> <.gff2> <.ortholog>\n";
	exit;
}

use strict;

my $gff1 = shift;
my $gff2 = shift;
my $orth = shift;

my (%chrGene1, %geneIndex);
&read_gff_ChrGene($gff1, \%chrGene1, \%geneIndex); # sub 1

my %chrGene2;
&read_gff_ChrGene($gff2, \%chrGene2, \%geneIndex); #sub 1

my %ortholog_geneSyn;
&gene_synteny_add($orth,\%ortholog_geneSyn,\%geneIndex,\%chrGene1); # sub 2

foreach my $chr(sort keys %ortholog_geneSyn)
{
	my @tmp = sort {$a->[2] <=> $b->[2]} @{$ortholog_geneSyn{$chr}};
	for (my $i = 0; $i < scalar @tmp; $i++)
	{
		my $ln;
		for (my $j = 0; $j < scalar @{$tmp[$i]};$j++)
		{
			$ln .= "$tmp[$i][$j]\t" 
		}
		$ln .= "\n";
		print $ln;
	}
}




## sub 1
##
sub read_gff_ChrGene
{

	my $gff_file = shift;
	my $chr_hshp = shift;
	my $ind_hshp = shift;

	open IN, $gff_file or die "Fail to open $gff_file\n";
	while (<IN>) 
	{
        	chomp;
        	my @info = split /\t/;
        	next unless ($info[2] eq 'mRNA');
        	die unless $info[8] =~ /ID=(\S+?);/;
        	my $gene = $1;
        	push @{$chr_hshp->{$info[0]}}, [$gene, $info[3], $info[4], $info[6]];
	}
	close IN;

	foreach my $chr (sort keys %{$chr_hshp}) 
	{
        	@{$chr_hshp->{$chr}} = sort {$a->[1] <=> $b->[1]} @{$chr_hshp->{$chr}};
        	for (my $i = 0; $i < @{$chr_hshp->{$chr}}; $i ++) 
		{
                	my ($gene, $bg, $ed, $strand) = @{$chr_hshp->{$chr}->[$i]};
                	my $index = $i + 1;
                	$ind_hshp->{$gene} = [$chr, $bg, $ed, $strand, $index];
        	}
	}

}

## sub 2
## $orth,\%ortholog_geneSyn,\%geneIndex,\%chrGene1
sub gene_synteny_add
{
	my $file = shift;
	my $hshp = shift;
	my $hsp2 = shift;
	my $hsp3 = shift;
	
	my %ind;
	open IN, $file or die "Fail to open $file\n";
	while (<IN>)
	{
		next if ($_ =~ m/#/);
		my ($gn1, $gn2, $level) = (split /\s+/)[0,1,-1];
		my ($g1,$g2) = ($gn1,$gn2);
		#$g1 =~ s/_[A-Z]{5}//;
		#$g2 =~ s/_[A-Z]{5}//;
		my ($chr1, $bg1, $ed1, $strand1, $index1) = @{$hsp2->{$g1}};
		die "$g2" unless ($hsp2->{$g2});
		my ($chr2, $bg2, $ed2, $strand2, $index2) = @{$hsp2->{$g2}};
		$ind{$chr1}{$index1} = 1;
		push @{$hshp->{$chr1}}, [$gn1,$chr1,$index1,$strand1,$bg1,$ed1,$gn2,$chr2,$index2,$strand2,$bg2,$ed2,$level];
	}
	close IN;

	foreach my $chr1 (keys %{$hsp3}) 
	{
	        foreach my $p (@{$hsp3->{$chr1}}) {
        	        my ($g1, $bg1, $ed1, $strand1) = @$p;
                	my $index1 = $hsp2->{$g1}->[-1];
	                next if ($ind{$chr1}{$index1});
        	        push @{$hshp->{$chr1}}, [$g1,$chr1,$index1,$strand1,$bg1,$ed1,"NA","NA","NA","NA","NA","NA","NA"];
        	}
	}
}


