#!/usr/bin/perl 
use strict;
use warnings;


die "Usage:$0 <genewise> <pep_len>\n\n" if @ARGV<1;

my $genewise=shift;
my $pep=shift;

my %Len;
open IN,$pep or die "$!";
while(<IN>){
	if (/(\S+)\s+(\S+)/){
		$Len{$1}=$2;
	}
}
close IN;

open IN,$genewise or die "$!";
$/="//\n";
while(<IN>){
	next if !/^Bits/;
	chomp;
	my $line1=$_;
	my $line2=<IN>;
	chomp $line2;
	my $line3=<IN>;	
	chomp $line3;
	my @l1=split(/\n/,$line1);
	my $len;
	my $pid;
	my $id;
	if (my @aa=split /\s+/,$l1[1]){
		$id=$aa[1];
		$pid=$id;
		$pid=$1 if($aa[1] =~ /(\S+)(-D\d+)/);
		$len=abs($aa[3]-$aa[2])+1;
	}else{
		die ;
	}
	my $cover=sprintf("%.2f",$len/$Len{$pid}*100);
	my @l3=split(/\n/,$line3);
	my $shift=-1;
	my ($frame,$memory);
	my @cds;
	my $chr;
	foreach (@l3){
		my @c=split(/\t/);
		my $start;
		if ($c[0]=~/^(\S+)_(\d+)_(\d+)/){
			$chr=$1;
			$start=$2;
		}else{
			die;
		}
		@c[3,4] =@c[4,3] if $c[3] >$c[4];
		if ($c[2] eq 'match'){
			$frame.="$memory\-".(($c[3]-1)+($start-1)).";" if($c[6] eq '+' && $shift >= 0);
			$frame.=(($c[4]+1)+($start-1))."-$memory".";"  if($c[6] eq '-' && $shift >= 0);
			$memory=($c[4]+1)+($start-1) if ($c[6] eq '+');
			$memory=($c[3]-1)+($start-1) if ($c[6] eq '-');
			$shift++;	
		}elsif ($c[2] eq 'cds'){
			$c[3]+=$start-1;
			$c[4]+=$start-1;
			push @cds,[$c[3],$c[4],$c[6],$c[7]];
		}
	}
	@cds=sort {$a->[0]<=>$b->[0]} @cds;
	my @mRNA=($frame)?($chr,'GeneWise','mRNA',$cds[0][0],$cds[-1][1],$cover,$cds[0][2],'.',"ID=$id;Shift=$shift;$frame"):($chr,'GeneWise','mRNA',$cds[0][0],$cds[-1][1],$cover,$cds[0][2],'.',"ID=$id;Shift=$shift;");
	print join("\t",@mRNA)."\n";
	foreach (@cds){
		print join("\t",$chr,'GeneWise','CDS',$$_[0],$$_[1],'.',$$_[2],$$_[3],"Parent=$id;")."\n";
	}
}
close IN;

