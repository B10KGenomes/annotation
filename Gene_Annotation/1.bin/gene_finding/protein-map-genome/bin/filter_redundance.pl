#!/usr/bin/perl
use strict;
use warnings;

#filter redundance of tab result(converted from net result).

die "Usage:$0 <tab> <outdir> <cutoff,default=100bp>\n\n" if @ARGV<2;

my $tab=shift;
my $Outdir=shift;
my $cutoff=shift;

$cutoff ||=100;
$Outdir=~s/\/$//;

my %Net;
read_tab($tab,\%Net);

my @Nr;
foreach my $id ( sort keys %Net ){
	@{$Net{$id}} = sort {$a->[7]<=>$b->[7]} @{$Net{$id}};
	foreach my $hsp ( @{$Net{$id}} ){
		if ( @Nr<1 ){
			push @Nr,[@$hsp];
		}elsif ( $Nr[-1][6] ne $id ){
			push @Nr,[@$hsp];
		}else{
			my $overlap=$Nr[-1][8]-$$hsp[7]+1;;
			if ( $overlap < $cutoff ){
				push @Nr,[@$hsp];
			}else{
				if ( $$hsp[9] > $Nr[-1][9] ){
					pop @Nr;
					push @Nr,[@$hsp];
				}else{
					next;
				}
			}
		}			
	}

}

@Nr = sort {$a->[0] cmp $b->[0] || $a->[2]<=>$b->[2]} @Nr;

my %Sort;
foreach my $hsp ( @Nr ){
	push @{$Sort{$hsp->[0]}},[@$hsp];
}

foreach my $chr ( sort keys %Sort ){
	open OUT,">$Outdir/$chr.nr.tab" or die "$!";
	foreach my $hsp ( @{$Sort{$chr}} ){
		print OUT join("\t",@$hsp)."\n";
	}
	close OUT;
}


# 0         1     2      3         4              5       6       7      8         9
#chr    chr_len start   end     align_len       strand  chr2    start2  end2    align_len2
sub read_tab{
	my ($file,$p)=@_;
	open IN,$file or die "$!";
	while(<IN>){
		next if /^#/;
		chomp;
		my @c=split(/\t/);
		push @{$p->{$c[6]}},[@c];	
	}
	close IN;
}
