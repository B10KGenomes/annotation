#!/usr/bin/perl

die "Usage:
	-s	extend number at start.[3]
	-e	extend number at end.[3]
	-t	threshold values.[10]

perl $0 [options] <gff> <out.id.lst> <genome.fa> > outfile\n\n" if @ARGV < 3;

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;

my ($add_s,$add_e,$over)=(3,3,10);
GetOptions(
	"s:i"=>\$add_s,
	"e:i"=>\$add_e,
	"t:i"=>\$over,
);

my $gfffile=shift;
my $lstfile=shift;
my $genomefile=shift;

#print STDERR "reading gff file ......\n";
my %gff;
read_gff($gfffile,\%gff);

#print STDERR "reading list file ......\n";
my %lst;
read_lst($lstfile,\%lst);


my %result;
open IN,$genomefile or die $!;
$/=">";<IN>;$/="\n";
while (<IN>){
	my $name=$1 if(/^(\S+)/);
#	print STDERR "dealing with scaffold $name ...\n";
	$/=">";
	chomp(my $seq=<IN>);
	$seq=~s/\s+//g;
	$/="\n";
	for my $id(keys %{$gff{$name}}){
#		print STDERR "dealing with gene $id ...\n";
		my $coor_s;
		my $coor_e;
		if($gff{$name}{$id}{mRNA}[0]->[6] eq '+'){
			$coor_s=$gff{$name}{$id}{mRNA}[0]->[3];
			$coor_e=$gff{$name}{$id}{mRNA}[0]->[4];
			for(my $i=0;$i <= $add_s+$lst{$id}->[0];$i++ ){
				last if($lst{$id}->[0] > $over);
				my $cor_s=$gff{$name}{$id}{mRNA}[0]->[3]-1;
				my $seq_s=substr($seq,$cor_s-$i*3,3);
				last if($seq_s =~ /n/i);
				next unless($seq_s =~ /ATG/i);
				$coor_s = $cor_s-$i*3+1;
				last;
			}
			for(my $i=0;$i <= $add_e+$lst{$id}->[1];$i++){
				last if($lst{$id}->[1] > $over);
				my $cor_e=$gff{$name}{$id}{mRNA}[0]->[4];
				my $seq_e=substr($seq,$cor_e+$i*3,3);
				last if($seq_e =~ /n/i);
				next unless($seq_e =~ /TAA/i || $seq_e =~ /TAG/i || $seq_e =~ /TGA/i);
				$coor_e = $cor_e+$i*3;
				last;
			}
		}else{
			$coor_s=$gff{$name}{$id}{mRNA}[0]->[4];
			$coor_e=$gff{$name}{$id}{mRNA}[0]->[3];
			for(my $i=0;$i <= $add_s+$lst{$id}->[0];$i++){
				last if($lst{$id}->[0] > $over);
				my $cor_s=$gff{$name}{$id}{mRNA}[0]->[4]-3;
				my $seq_s=substr($seq,$cor_s+$i*3,3);
				$seq_s=reverse $seq_s;
				$seq_s =~ tr/atgcATGC/tacgTACG/;
				last if($seq_s =~ /n/i);
				next unless($seq_s =~ /ATG/i);
				$coor_s = $cor_s+$i*3+3;
				last;
			}
			for(my $i=0;$i <= $add_e+$lst{$id}->[1];$i++){
				last if($lst{$id}->[1] > $over);
				my $cor_e=$gff{$name}{$id}{mRNA}[0]->[3]-4;
				my $seq_e=substr($seq,$cor_e-$i*3,3);
				$seq_e=reverse $seq_e;
				$seq_e =~ tr/atgcATGC/tacgTACG/;
				last if($seq_e =~ /n/i);
				next unless($seq_e =~ /TAA/i || $seq_e =~ /TAG/i || $seq_e =~ /TGA/i);
				$coor_e = $cor_e+4;
				last;
			}
		}
		@{$result{$id}}=($coor_s,$coor_e);
#		print STDERR "deal with gene $id over !\n";
	}
#	print STDERR "read scaffold $name over ! \n";
}
close IN;

#print Dumper %result;

for my $name(keys %gff){
	for my $id(keys %{$gff{$name}}){
		@{$gff{$name}{$id}{cds}}=sort {$a->[3] <=> $b->[4]} @{$gff{$name}{$id}{cds}};
		my ($start,$end)=(0,0);
		if($result{$id}->[0] <= $result{$id}->[1]){
			$start=1 if($gff{$name}{$id}{cds}[0]->[3] != $result{$id}->[0]);
			$end=1 if($gff{$name}{$id}{cds}[-1]->[4] != $result{$id}->[1]);
			$gff{$name}{$id}{cds}[0]->[3]=$result{$id}->[0];
			$gff{$name}{$id}{cds}[-1]->[4]=$result{$id}->[1];
		}else{
			$start=1 if($gff{$name}{$id}{cds}[-1]->[4] != $result{$id}->[0]);
			$end=1 if($gff{$name}{$id}{cds}[0]->[3] != $result{$id}->[1]);
			$gff{$name}{$id}{cds}[0]->[3]=$result{$id}->[1];
			$gff{$name}{$id}{cds}[-1]->[4]=$result{$id}->[0];
		}
		($result{$id}->[0],$result{$id}->[1])=($result{$id}->[1],$result{$id}->[0]) if($result{$id}->[0] > $result{$id}->[1]);
		print STDERR "$id\t$start\t$end\t$gff{$name}{$id}{mRNA}[0]->[6]\t$gff{$name}{$id}{mRNA}[0]->[3]\-$result{$id}->[0]\t$gff{$name}{$id}{mRNA}[0]->[4]\-$result{$id}->[1]\n" if($start ==1 || $end ==1);
		($gff{$name}{$id}{mRNA}[0]->[3],$gff{$name}{$id}{mRNA}[0]->[4])=($result{$id}->[0],$result{$id}->[1]);
		print join("\t",@{$gff{$name}{$id}{mRNA}[0]})."\n";
		for(@{$gff{$name}{$id}{cds}}){
			print join("\t",@{$_})."\n";
		}
	}
}


######################
sub read_gff{
	my ($file,$hash)=@_;
	open IN,$file or die $!;
	while (<IN>){
		next if /^#/;
		chomp;
		my @a=split /\t+/;
		my $id=$2 if ($a[8]=~/(ID=|Parent=)([^;\s]+)/);
#		print STDERR "read in $id ...\n";
		($a[3],$a[4])=($a[4],$a[3]) if($a[3] > $a[4]);
		if($a[2]=~/mRNA/){
			push @{$$hash{$a[0]}{$id}{mRNA}},[@a];
		}else{
			push @{$$hash{$a[0]}{$id}{cds}},[@a];
		}
	}
	close IN;
}

sub read_lst{
	my ($file,$hash)=@_;
	open IN,$file or die $!;
	while (<IN>){
		next if /^#/;
		chomp;
		my @a=split /\t+/;
		@{$$hash{$a[0]}}=($a[8],$a[9]);
	}
	close IN;
}
