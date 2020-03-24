#!/usr/bin/perl

use strict;

die "perl $0 <gff> > outfile \n" if @ARGV < 1;

my $Ref=shift;
my $File=$Ref;

my %gff;
read_gff($File,\%gff);

my %ref;
read_gff($Ref,\%ref);

for my $name(keys %ref){
	for my $id(sort {$ref{$name}{$a}{mRNA}->[3] <=> $ref{$name}{$b}{mRNA}->[3]} keys %{$ref{$name}}){
		my $q=$ref{$name}{$id};
		my @aa=sort {$gff{$name}{$a}{mRNA}->[3] <=> $gff{$name}{$b}{mRNA}->[3]} keys %{$gff{$name}};
		for(@aa){
			my $p=$gff{$name}{$_};
			next if($p->{mRNA}->[4] < $q->{mRNA}->[3]);
			last if($p->{mRNA}->[3] > $q->{mRNA}->[4]);
			my $mRNA_over=overlap($q->{mRNA}->[3],$q->{mRNA}->[4],$p->{mRNA}->[3],$p->{mRNA}->[4]);
			my @q_cds=sort {$a->[3] <=> $b->[3]} @{$ref{$name}{$id}{CDS}};
			my @p_cds=sort {$a->[3] <=> $b->[3]} @{$gff{$name}{$_}{CDS}};
			my ($len_q_cds,$len_p_cds,$cds_over)=(0,0,0);
			for(@q_cds){
				$len_q_cds+=($_->[4]-$_->[3]+1);
			}
			for(@p_cds){
				$len_p_cds+=($_->[4]-$_->[3]+1);
			}
			for my $i(@q_cds){
				for my $j(@p_cds){
					next if($j->[4] < $i->[3]);
					last if($j->[3] > $i->[4]);
					$cds_over+=overlap($i->[3],$i->[4],$j->[3],$j->[4]);
				}
			}
			my $short = ($len_q_cds < $len_p_cds) ? $len_q_cds : $len_p_cds;
			my $cds_dist = 1 - $cds_over/$short;
			print "$id\t$_";
			printf "\t%.4f\t%g", $cds_dist,$cds_over;
			print "\t$name";
			print "\t$ref{$name}{$id}{mRNA}->[6]\t$len_q_cds\t$ref{$name}{$id}{mRNA}->[3]\t$ref{$name}{$id}{mRNA}->[4]";
			print "\t$gff{$name}{$_}{mRNA}->[6]\t$len_p_cds\t$gff{$name}{$_}{mRNA}->[3]\t$gff{$name}{$_}{mRNA}->[4]\n";
		}
	}
}


############
sub overlap{
	my ($aa,$ab,$ba,$bb)=@_;
	my $max=$ab>$bb?$ab:$bb;
	my $min=$aa<$ba?$aa:$ba;
	my $ov=0;
	$ov=($bb-$ba+$ab-$aa+2)-($max-$min+1);
	$ov=0 if ($ov < 0);
	return $ov;
}
sub read_gff{
	my ($file,$hash)=@_;
	open IN,$file or die $!;
	while(<IN>){
		next if /^#/;
		chomp;
		my @a=split /\t+/;
		my $name=$a[0];
		@a[3,4]=@a[4,3] if ($a[3] > $a[4]);
		my $id;
		if($a[2] eq 'mRNA' && $a[8] =~ /ID=([^;\s]+)/){
			$id=$1;
			@{$$hash{$name}{$id}{mRNA}}=@a;
		}elsif(($a[2] eq 'CDS' || $a[2] eq 'exon')&& $a[8] =~ /Parent=([^;\s]+)/){
			$id=$1;
			push @{$$hash{$name}{$id}{CDS}},[@a];
		}
	}
	close IN;
}
