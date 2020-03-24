#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($Bin $Script);

die "Usage:$0 <id.lst> <pep1.pep (advise input genewise result)> <pep2.pep>\n" if @ARGV<2;

my $lst = shift;
my $fa1 = shift;
my $fa2 = shift;

my %Pep1;
Read_fasta($fa1,\%Pep1);

my %Pep2;
Read_fasta($fa2,\%Pep2);

my %config;
parse_config("$Bin/config.txt",\%config);

my $sh_out = '';
#my $name = `basename $lst`;
my $name = $lst;
chomp $name;
my $pwd = $ENV{'PWD'};
my $dir = "$name.cut".time();
`rm -r $dir` while -e $dir;
`mkdir $dir` unless -d $dir;

my $loop2 = 0;
my $loop1 = 0;
my $subdir2 = "000";
my $subdir1 = "000";

open IN,$lst or die "$!";
open OUTSH,">$name.muscle.sh" or die "$!";
while(<IN>){
	my ($id1,$id2);
	if (/^(\S+)\s+(\S+)/){
		$id1 = $1;
		$id2 = $2;
	}else{
		die "List format error in line $.\n";
	}
	if ( (not defined $id1) or (not defined $id2) ){
		die "Protein sequence of $id1/$id2 is not found\n";
	}
	if($loop2 % 200 == 0){
		if($loop1 % 100 == 0){
			$subdir1++;
			mkdir ("$dir/$subdir1") unless -d "$dir/$subdir1";
			$subdir2 = "000";
		 }
		 $subdir2++;
		 mkdir("$dir/$subdir1/$subdir2") unless -d "$dir/$subdir1/$subdir2";
		 $loop1++;
	}
	my $sub_pep = "$dir/$subdir1/$subdir2/$id1.$id2.fa";
	my $muscle = "$dir/$subdir1/$subdir2/$id1.$id2.muscle";
	my $log="$dir/$subdir1/$subdir2/$id1.$id2.log";

	open OUT,">$sub_pep" or die "$!";
	print OUT ">$id1\n$Pep1{$id1}\n".">$id2\n$Pep2{$id2}\n";
	close OUT;

	$sh_out = "$config{muscle} -in $sub_pep -out $muscle 2> $log;\n";
	print OUTSH $sh_out;
	$loop2++;
}
close IN;

#open OUT,">$name.muscle.sh" or die "$!";
#print OUT $sh_out;
close OUTSH;


#read fasta file
#usage: Read_fasta($file,\%hash);
#############################################
sub Read_fasta{
	my $file=shift;
	my $hash_p=shift;
	
	my $total_num;
	open(IN, $file) || die ("can not open $file\n");
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		chomp;
		my ($name,$seq);
		
		if ( /^(\S+)/ ){
			$name = $1;
		}else{
			die "Fasta format error!\n";
		}

		$/=">";
		$seq = <IN>;
		chomp $seq;
		$seq=~s/\s//g;
		$/="\n";
		
		$hash_p->{$name} = $seq;

		$total_num++;
	}
	close(IN);
	
	return $total_num;
}

sub parse_config{
    my $conifg_file = shift;
    my $config_p = shift;

    my $error_status = 0;
    open IN,$conifg_file || die "fail open: $conifg_file";
    while (<IN>) {
        next if /^#/;
        if (/(\S+)\s*=\s*(\S+)/) {
            my ($software_name,$software_address) = ($1,$2);
            $config_p->{$software_name} = $software_address;
            if (! -e $software_address){
                warn "Non-exist:  $software_name  $software_address\nPlease re-configure the config.txt file\n";
                $error_status = 1;
            }
        }
    }
    close IN;
    die "\nExit due to error of software configuration\n" if($error_status);
}
