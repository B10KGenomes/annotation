#!/usr/bin/perl -w

=head1 Description

This script is used to run the RBH, but only the L1 level is the true RBH, so please grep L1 from th results.

=head1 Contact & Version

  Author: Yuan Deng & Qi Fang, dengyuan@genomics.cn & fangqi@genomics.cn
  Version: p1,   Date: 2019-11-19
  Version: p2,   Date: 2020-01-17 Modified by Qi Fang to choose the best synteny

=head1 Usage Exmples

  perl 97.ort_rbh.syn.pl intersect.ortholog.list.order.ort.cor reference_gene_tag query_gene_tag > results

=cut

use strict;
unless (@ARGV == 3) {
	die <<End;

Description:
    This script is used to get ortholog pairs.

Usage:
    perl $0 <idAdd> <ref_species> <hit_species> > ref_hit.idAdd.ort

End
}
my $tab_file = shift;
my $species1 = shift;
my $species2 = shift;
my @speciesOrder = ($species1, $species2);
##################################################################################################################

## get the gene pairs which are best hit to each other.
my %hit; ## mark the hits of a query
my %bestPair; ## mark the RBH pairs
&getRBH($tab_file, $species1, $species2, \%hit, \%bestPair);

my @result; ## rank the hit pairs
my %genePair; ## mark the gene pairs determine by mcscan.
## output result and mark the gene pairs added.
foreach my $g1 (keys %bestPair) {
	foreach my $g2 (keys %{$bestPair{$g1}}) {
        next if $bestPair{$g1}{$g2}[3] <= 50;
		next if $genePair{$g1} || $genePair{$g2};
		($g1, $g2) = &getGeneOrder($g1, $g2);
		my ($align1, $align2, $align, $score) = @{$bestPair{$g1}{$g2}};
		push @result, "$g2\t$g1\tL1\t$score\n";
		$genePair{$g1} ++;
		$genePair{$g2} ++;
	}
}

## output result

#print "#Referrence \tSubject \tR_ratio \tS_ratio \tRatio \tSocre \tLevel \n";  
print "$_" foreach @result;


##subroutine
##################################################################
sub getRBH {
my ($tab_file, $species1, $species2, $hit_ref, $bestPair_ref) = @_;
my %species = ($species1=>1, $species2=>1);
## restore the hits for each query
my %hit;
open IN, $tab_file;
while (<IN>) {
	next if /^#/;
	my ($g1, $align1, $g2, $align2, $score) = (split /\s+/)[0,4,6,10,11];
	next if $g1 eq $g2;
	my $sp1 = (split /_/, $g1)[0];
	my $sp2 = (split /_/, $g2)[0];
	next unless $species{$sp1} && $species{$sp2};
	my $align = ($align1 + $align2) / 2; ## take the average alignRate for sort
	push @{$hit{$g1}}, [$g2, $align1, $align2, $align, $score];
	push @{$hit_ref->{$g1}}, [$g2, $align1, $align2, $align, $score];
	push @{$hit_ref->{$g2}}, [$g1, $align2, $align1, $align, $score];
}
close IN;
## mark the best hit for each query
my %pair;
my %links; ##designed for ALL
foreach my $g1 (keys %hit) {
	@{$hit{$g1}} = sort {$b->[-1] <=> $a->[-1]} @{$hit{$g1}};
##filter out any group with Tandem
    if($hit{$g1}[-1][4] < 50 or ($hit{$g1}[1] and $hit{$g1}[1][4]==$hit{$g1}[0][4]))
    {##if singleton connecting Perfect would better be single copy orthologs But them are rare
        next, unless($g1=~/$species1/);
        if(%links){
            my $overlaped=0;
            foreach my $tt (sort keys %links) {
                foreach my $p (@{$hit{$g1}}) {
                    if(grep { $_ eq ${$p}[0] } @{$links{$tt}}){
                        $overlaped++;
                        last;
                    }
                }

                if($overlaped>0){
                    my $ttt = $tt.",".$g1;
                    @{$links{$ttt}}=@{$links{$tt}};
                    delete $links{$tt};
                    foreach my $p (@{$hit{$g1}}) {
                        push @{$links{$ttt}}, ${$p}[0];
                    }
                    last;
                }
            }

            if($overlaped==0){
                foreach my $p (@{$hit{$g1}}) {
                    push @{$links{$g1}}, ${$p}[0];
                }
            }

        }else{
            foreach my $p (@{$hit{$g1}}) {
                push @{$links{$g1}}, ${$p}[0];
            }
        }
    }
    else
    {
	    foreach my $p (@{$hit{$g1}}) {
#            if($hit{$g1}[1] and $hit{$g1}[1][4]==$hit{$g1}[0][4])
#            {##segement duplications
#                print "$g1\t$hit{$g1}[0][0]:$hit{$g1}[0][4]\t$hit{$g1}[1][0]:$hit{$g1}[1][4]\n";
#            }
#            else
#            {
        		my ($g2, $align1, $align2, $align, $score) = @$p;
	        	$pair{$g1}{$g2} = [$align1, $align2, $align, $score];
#            }
		    last;
    	}
    }
}
## restore RBH pairs
my %printMark;
foreach my $g1 (sort keys %pair) {
	foreach my $g2 (keys %{$pair{$g1}}) {
		next unless $pair{$g2}{$g1}; ##make sure both best
		($g1, $g2) = sort ($g1, $g2);

		my ($align1_1, $align1_2, $align11, $score11) = @{$pair{$g1}{$g2}};
		my ($align2_1, $align2_2, $align22, $score22) = @{$pair{$g2}{$g1}};
		my ($align1, $align2, $align, $score);
		if ($score22 > $score11) {
			($align1, $align2, $align, $score) = ($align2_2, $align2_1, $align22, $score22);	
		} else {
			($align1, $align2, $align, $score) = ($align1_1, $align1_2, $align11, $score11);
		}
		unless ($printMark{$g1}{$g2}) {
			$bestPair_ref->{$g1}{$g2} = [$align1, $align2, $align, $score];
			$bestPair_ref->{$g2}{$g1} = [$align2, $align1, $align, $score];
		}
		$printMark{$g1}{$g2} ++;
	}
}

## print non-synteny orthologous groups
foreach my $tt (sort keys %links) {
    my @a = split(/,/,$tt);
    my @g1 = sort {$a cmp $b} @a;
    my %temp;
    my @temp=grep {++$temp{$_}<2} @{$links{$tt}};
    my @g2=sort {$a cmp $b} @temp;
    print join(",",@g2)."\t".join(",",@g1)."\tL9\t40\n";
}
}

###############################
sub getGeneOrder {
my ($g1, $g2) = @_;
my $sp1 = (split /_/, $g1)[0];
my $sp2 = (split /_/, $g2)[0];
my %sp_to_gene = ($sp1=>$g1, $sp2=>$g2);
($g1, $g2) = ($sp_to_gene{$speciesOrder[0]}, $sp_to_gene{$speciesOrder[1]});
return ($g1, $g2);
}
