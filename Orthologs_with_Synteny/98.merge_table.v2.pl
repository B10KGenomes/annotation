#!/usr/bin/perl -w

=head1 Description

This script is used to merge all the orthologs to one table.

=head1 Contact & Version

  Author: Yuan Deng & Qi Fang, dengyuan@genomics.cn & fangqi@genomics.cn
  Version: p1,   Date: 2019-11-19

=head1 Usage Exmples

  perl 98.merge_table.pl input_list_file results

=cut


if (@ARGV != 2)
{
	print "perl $0 <ortholog_files.list> <prefix for output: ortholog>\n";
	exit;
}

use strict;
use File::Basename qw(basename dirname);

my $list = shift;
my $out  = shift;
my $abb_file = "/hwfssz5/ST_DIVERSITY/B10K/USER/dengyuan/01.assembly/short_name";

my %ortholog;
&read_files($list,\%ortholog); # sub 1

my @taxa;
&taxa_from_files($list,\@taxa); # sub 2

my %short;
&short_name($abb_file,\%short); # sub 3

my @taxa2 = @taxa;
shift @taxa2;

##merge the ref of %ortholog
#get the ref list
my @tmpKeys = sort keys %ortholog;
my %pairKeys;
foreach my $k (@tmpKeys){
    my @a = split(/,/,$k);
    for(my$i=0;$i<=$#a;$i++){
        push @{$pairKeys{$a[$i]}},@a;
    }
}
#filter ref list
my %pairKeys2;
foreach my $yyu (sort keys %pairKeys){
    my %temp;
    my @temp1=grep {++$temp{$_}<2} @{$pairKeys{$yyu}};
    my @temp2=sort {$a cmp $b} @temp1;
    $pairKeys{$yyu}=\@temp2;
}
#re-merge the ref list
my @tmpHash;
foreach my $yyu (sort keys %pairKeys){
    foreach my $yyy (sort {$a cmp $b} @{$pairKeys{$yyu}}){
        push @{$pairKeys2{$yyu}}, @{$pairKeys{$yyy}};
    }
    my %temp;
    my @temp1=grep {++$temp{$_}<2} @{$pairKeys2{$yyu}};
    my @temp2=sort {$a cmp $b} @temp1;
    $pairKeys2{$yyu}=\@temp2;
    push @tmpHash, join(",",@temp2);
}
#re-filter the ref list
my %tmpHashtmp;
my @tempHashtmp1=grep {++$tmpHashtmp{$_}<2} @tmpHash;
@tmpHash=@tempHashtmp1;
print join("\n",@tmpHash)."\n";
#filter %ortholog based on new ref list
my %ortholog2;
foreach my $ref (sort keys %ortholog) {
    my @a = split(/,/,$ref);
    foreach my $ref2 (@tmpHash){
        my @b = split(/,/,$ref2);
        if(grep {$_ eq $a[0]} @b){
            foreach (@taxa2){
                my $abb_name=$short{$_};
                if (exists $ortholog{$ref}{$abb_name} and exists $ortholog2{$ref2}{$abb_name}){
                    $ortholog2{$ref2}{$abb_name}.=",".$ortholog{$ref}{$abb_name};
                }elsif(exists $ortholog{$ref}{$abb_name}){
                    $ortholog2{$ref2}{$abb_name}=$ortholog{$ref}{$abb_name};
                }
            }
        }
    }
}

my $i = 1;
my ($out1,$out2);
foreach my $ref(sort keys %ortholog2)
{
	$out1 .= "$i\t$ref\t";
	$out2 .= "$i\t1\t";
	my $num = 1;
	foreach (@taxa2)
	{
		my $abb_name=$short{$_};
		if (exists $ortholog2{$ref}{$abb_name})
		{
			$out2 .= "1\t";
			$out1 .= "$ortholog2{$ref}{$abb_name}\t";
			$num ++;
		}
		else
		{
			$out2 .= "0\t";
			$out1 .= "-\t"
		}
	}
	$out1 .= "\n";
	$out2 .= "$num\n";
	$i++;
}

open OUT, ">${out}_ortholog_group" or die;
print OUT "ID\t".join ("\t", @taxa)."\n";
print OUT $out1;
close OUT;

open OUT, ">${out}_ortholog_group.stat" or die;
print OUT "ID\t".join ("\t", @taxa)."\tGeneNum\n";
print OUT $out2;
close OUT;




## sub 1
##
sub read_files
{
	my $file = shift;
	my $hshp = shift;

	open IN,$file or die "Fail to open $file\n";
	while (my $ln = <IN>)
	{
		chomp $ln;
		open ORTH,$ln or die "Fail to open $ln\n";
		while (<ORTH>)
		{
			my ($id1,$id2) = (split)[1,0];
			$id2 =~ /^([\w-]+)_/;
			$hshp->{$id1}{$1} = $id2;
		}
		close ORTH;
		
	}
	close IN;
}


## sub 2
##
sub taxa_from_files
{
	my $file = shift;
	my $hshp = shift;

	open IN,$file or die;
	my $head = <IN>;
	my $taxa = basename($head);
	$taxa =~ /([A-Z]{6})_([\w,-]+).*.ort/;
	push @$hshp,$1;
	push @$hshp,$2;
	while (<IN>)
	{
		$taxa = basename($_);
		$taxa =~ /[A-Z]{6}_([\w,-]+).*.ort/;
		push @$hshp,$1;
	}
	close IN;
}

## sub 3
##
sub short_name
{
	my $file = shift;
	my $hshp = shift;
	open IN,$file or die;
	while (<IN>){
		chomp;
		$_ =~ /([\w,-]+)\s+([A-Z]{6})/;
		$hshp->{$1}=$2;
	}
	close IN;
}
