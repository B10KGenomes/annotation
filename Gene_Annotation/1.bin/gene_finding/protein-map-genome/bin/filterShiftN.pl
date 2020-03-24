#!/usr/bin/perl

die "perl $0 <gff> <genome.fa> [extend-extend] > outfile\n" if(@ARGV < 2);

use strict;

my $gfffile=shift;
my $genomefile=shift;
my $extend=shift;
$extend ||= "6-6";
my ($ex_s,$ex_e)=($1,$2) if ($extend =~ /^(\d+)-(\d+)/);

my %gff;
read_gff($gfffile,\%gff);
#my %fa;
#read_fa($genomefile,\%fa);


open IN,$genomefile or die $!;
$/=">";<IN>;$/="\n";
while(<IN>){
	my $name=$1 if(/^(\S+)/);
	$/=">";
	chomp(my $seqq=<IN>);
	$/="\n";
	$seqq=~s/\s+//g;
	for my $id(keys %{$gff{$name}}){
		my $shift=$gff{$name}{$id};
		my @a=split /;/,$shift;
		my $flag=0;
		for(@a){
			my ($s,$e)=($1,$2) if /(\d+)-(\d+)/;
			my $seq=substr($seqq,$s-1-$ex_s,$e-$s+1+$ex_s+$ex_e);
			$flag++ if($seq !~ /N/i);
		}
		print "$id\n" if($flag == 0);
	}
}

################
sub read_gff{
	my ($file,$hash)=@_;
	open IN,$file or die $!;
	while (<IN>){
		next if /^#/;
		next unless /mRNA/;
		next if /Shift=0/i;
		chomp;
		my @a=split /\t+/;
		my $sign=$a[6];
		if($a[8]=~/ID=(\S+);Shift=(\d+);(\S+);/){
			my $shift=$3;
			my $id=$1;
			$$hash{$a[0]}{$id}=$shift;
		}
	}
	close IN;
}

sub read_fa{
	my ($file,$hash)=@_;
	open IN,$file or die $!;
	$/=">";<IN>;$/="\n";
	while(<IN>){
		my $id=$1 if(/^(\S+)/);
		$/=">";
		chomp(my $seq=<IN>);
		$/="\n";
		$seq=~s/\s+//g;
		$$hash{$id}=$seq;
	}
	close IN;
}
