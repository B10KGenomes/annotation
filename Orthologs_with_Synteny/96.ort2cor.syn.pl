#!/usr/bin/perl -w

=head1 Description

This script is used to translate the synteny information to follower steps.

=head1 Contact & Version

  Author: Qi Fang, fangqi@genomics.cn
  Version: p1,   Date: 2020-01-17

=head1 Usage Exmples

  perl 96.ort2cor.syn.pl cyc2tandem intersect.ortholog.list.order.ort > intersect.ortholog.list.order.ort.cor

=cut

use strict;
my $cyc = shift;
my $input = shift;

my (%chra,%chrb); ##store the order
my (%genea,%geneb); ##store the gene
my %pairs; ##store the synteny relation
my (%hasha,%hashb);
open IN,$input;
while(<IN>){
    chomp;
    my @a = split /\t+/;
    push @{$chra{$a[1]}}, $a[2];
    push @{$chrb{$a[7]}}, $a[8];
    $genea{$a[1]}{$a[2]}=$a[0];
    $geneb{$a[7]}{$a[8]}=$a[6];
    $pairs{$a[0]}{$a[6]}=$a[3]."\t".$a[9]."\t".$a[12]."\t".$a[13];
    $pairs{$a[6]}{$a[0]}=$a[9]."\t".$a[3]."\t".$a[12]."\t".$a[13];
    if(!exists $hasha{$a[0]}){
        push @{$hasha{$a[0]}}, ($a[1],$a[2]);
    }
    if(!exists $hashb{$a[6]}){
        push @{$hashb{$a[6]}}, ($a[7],$a[8]);
    }
}close IN;

##sorting and searching
my (%tagSa,%tagSb); ##single copy
my (%tagDa,%tagDb); ##tendam copy

foreach my $k (sort keys %chra){
    my %temp;
    my @temp1=grep {++$temp{$_}<2} @{$chra{$k}};
    my @temp2=sort {$a <=> $b} @temp1;
    $chra{$k}=\@temp2;

    foreach my $g (@{$chra{$k}}){
        if(exists $pairs{$genea{$k}{$g}}){
            my $len = %{$pairs{$genea{$k}{$g}}};
            if($len==1){
                $tagSa{$genea{$k}{$g}}=1;$tagDa{$genea{$k}{$g}}=0;next;
            }elsif($len<1){
                print "$genea{$k}{$g}\tError\n";
            }else{
                $tagSa{$genea{$k}{$g}}=0;
                for(my $i=1;$i<=$cyc;$i++){
###cycle for tendam length
                    my $last;
                    foreach my $t (sort{$hashb{$a}[0] cmp $hashb{$b}[0] or $hashb{$a}[1] <=> $hashb{$b}[1]} keys %{$pairs{$genea{$k}{$g}}}){
                        last, if($i>1 and $tagDa{$genea{$k}{$g}}>=1);#####TEST
                        unless($last){$last=$hashb{$t}[0];$tagDa{$genea{$k}{$g}}=0;next;}
                        if($last eq $hashb{$t}[0] and ($hashb{$t}[1]-$i)>0 and exists $geneb{$last}{$hashb{$t}[1]-$i} and exists $pairs{$genea{$k}{$g}}{$geneb{$last}{$hashb{$t}[1]-$i}})
                        {
                            if($i>1){
                                my $tmp_multiDup=0;
                                my $str;
                                for(my $j=1; $j<=$i-1; $j++){
                                    if($j==1){
                                        if($g-$j>0 and ($hashb{$t}[1]-$i-$j)>0 and exists $genea{$k}{$g-$j} and exists $geneb{$last}{$hashb{$t}[1]-$j} and exists $geneb{$last}{$hashb{$t}[1]-$i-$j} and exists $pairs{$genea{$k}{$g-$j}}{$geneb{$last}{$hashb{$t}[1]-$j}} and exists $pairs{$genea{$k}{$g-$j}}{$geneb{$last}{$hashb{$t}[1]-$i-$j}})
                                        {
                                            $tmp_multiDup=1; $str=1;
                                        }
                                        elsif($g-$j>0 and ($hashb{$t}[1]-$i-$j)>0 and exists $genea{$k}{$g+$j} and exists $geneb{$last}{$hashb{$t}[1]-$j} and exists $geneb{$last}{$hashb{$t}[1]-$i-$j} and exists $pairs{$genea{$k}{$g+$j}}{$geneb{$last}{$hashb{$t}[1]-$j}} and exists $pairs{$genea{$k}{$g+$j}}{$geneb{$last}{$hashb{$t}[1]-$i-$j}})
                                        {
                                            $tmp_multiDup=1; $str=2;
                                        }
                                        else{last;}
                                    }else{
                                        if($str==1){
                                            unless($g-$j>0 and ($hashb{$t}[1]-$i-$j)>0 and exists $genea{$k}{$g-$j} and exists $geneb{$last}{$hashb{$t}[1]-$j} and exists $geneb{$last}{$hashb{$t}[1]-$i-$j} and exists $pairs{$genea{$k}{$g-$j}}{$geneb{$last}{$hashb{$t}[1]-$j}} and exists $pairs{$genea{$k}{$g-$j}}{$geneb{$last}{$hashb{$t}[1]-$i-$j}})
                                            {$tmp_multiDup=0;}
                                        }elsif($str==2){
                                            unless($g-$j>0 and ($hashb{$t}[1]-$i-$j)>0 and exists $genea{$k}{$g+$j} and exists $geneb{$last}{$hashb{$t}[1]-$j} and exists $geneb{$last}{$hashb{$t}[1]-$i-$j} and exists $pairs{$genea{$k}{$g+$j}}{$geneb{$last}{$hashb{$t}[1]-$j}} and exists $pairs{$genea{$k}{$g+$j}}{$geneb{$last}{$hashb{$t}[1]-$i-$j}})
                                            {$tmp_multiDup=0;}
                                        }
                                    }
                                }
                                if($tmp_multiDup==1){$tagDa{$genea{$k}{$g}}=$i;}
                            }else{
                                $tagDa{$genea{$k}{$g}}=1;
                            }
                        }elsif($last eq $hashb{$t}[0] and exists $geneb{$last}{$hashb{$t}[1]+$i} and exists $pairs{$genea{$k}{$g}}{$geneb{$last}{$hashb{$t}[1]+$i}})
                        {
                            if($i>1){
                                my $tmp_multiDup=0;
                                my $str;
                                for(my $j=1; $j<=$i-1; $j++){
                                    if($j==1){
                                        if(exists $genea{$k}{$g+$j} and exists $geneb{$last}{$hashb{$t}[1]+$j} and exists $geneb{$last}{$hashb{$t}[1]+$i+$j} and exists $pairs{$genea{$k}{$g+$j}}{$geneb{$last}{$hashb{$t}[1]+$j}} and exists $pairs{$genea{$k}{$g+$j}}{$geneb{$last}{$hashb{$t}[1]+$i+$j}})
                                        {
                                            $tmp_multiDup=1;$str=1;
                                        }
                                        elsif($g-$j>0 and exists $genea{$k}{$g-$j} and exists $geneb{$last}{$hashb{$t}[1]+$j} and exists $geneb{$last}{$hashb{$t}[1]+$i+$j} and exists $pairs{$genea{$k}{$g-$j}}{$geneb{$last}{$hashb{$t}[1]+$j}} and exists $pairs{$genea{$k}{$g-$j}}{$geneb{$last}{$hashb{$t}[1]+$i+$j}})
                                        {
                                            $tmp_multiDup=1;$str=2;
                                        }else{last;}
                                    }else{
                                        if($str==1){
                                            unless(exists $genea{$k}{$g+$j} and exists $geneb{$last}{$hashb{$t}[1]+$j} and exists $geneb{$last}{$hashb{$t}[1]+$i+$j} and exists $pairs{$genea{$k}{$g+$j}}{$geneb{$last}{$hashb{$t}[1]+$j}} and exists $pairs{$genea{$k}{$g+$j}}{$geneb{$last}{$hashb{$t}[1]+$i+$j}})
                                            {$tmp_multiDup=0;}
                                        }elsif($str==2){
                                            unless($g-$j>0 and exists $genea{$k}{$g-$j} and exists $geneb{$last}{$hashb{$t}[1]+$j} and exists $geneb{$last}{$hashb{$t}[1]+$i+$j} and exists $pairs{$genea{$k}{$g-$j}}{$geneb{$last}{$hashb{$t}[1]+$j}} and exists $pairs{$genea{$k}{$g-$j}}{$geneb{$last}{$hashb{$t}[1]+$i+$j}})
                                            {$tmp_multiDup=0;}
                                        }
                                    }
                                }
                                if($tmp_multiDup==1){$tagDa{$genea{$k}{$g}}=$i;}
                            }else{
                                $tagDa{$genea{$k}{$g}}=1;
                            }
                        }else{
                            $tagDa{$genea{$k}{$g}}=0,unless($tagDa{$genea{$k}{$g}}>=1);
                            $last=$hashb{$t}[0];
                        }
                    }
###cycle for tendam length
                }
            }
        }
    }
}

foreach my $k (sort keys %chrb){
    my %temp;
    my @temp1=grep {++$temp{$_}<2} @{$chrb{$k}};
    my @temp2=sort {$a <=> $b} @temp1;
    $chrb{$k}=\@temp2;

    foreach my $g (@{$chrb{$k}}){
        if(exists $pairs{$geneb{$k}{$g}}){
            my $len = %{$pairs{$geneb{$k}{$g}}};
            if($len==1){
                $tagSb{$geneb{$k}{$g}}=1;$tagDb{$geneb{$k}{$g}}=0;next;
            }elsif($len<1){
                print "$geneb{$k}{$g}\tError\n";
            }else{
                $tagSb{$geneb{$k}{$g}}=0;
                for(my $i=1;$i<=$cyc;$i++){
###cycle for tendam length
                    my $last;
                    foreach my $t (sort{$hasha{$a}[0] cmp $hasha{$b}[0] or $hasha{$a}[1] <=> $hasha{$b}[1]} keys %{$pairs{$geneb{$k}{$g}}}){
                        last, if($i>1 and $tagDb{$geneb{$k}{$g}}>=1);#####TEST
                        unless($last){$last=$hasha{$t}[0];$tagDb{$geneb{$k}{$g}}=0;next;}
                        if($last eq $hasha{$t}[0] and ($hasha{$t}[1]-$i)>0 and exists $genea{$last}{$hasha{$t}[1]-$i} and exists $pairs{$geneb{$k}{$g}}{$genea{$last}{$hasha{$t}[1]-$i}})
                        {
                            if($i>1){
                                my $tmp_multiDup=0;
                                my $str;
                                for(my $j=1; $j<=$i-1; $j++){
                                    if($j==1){
                                        if($g-$j>0 and ($hasha{$t}[1]-$i-$j)>0 and exists $geneb{$k}{$g-$j} and exists $genea{$last}{$hasha{$t}[1]-$j} and exists $genea{$last}{$hasha{$t}[1]-$i-$j} and exists $pairs{$geneb{$k}{$g-$j}}{$genea{$last}{$hasha{$t}[1]-$j}} and exists $pairs{$geneb{$k}{$g-$j}}{$genea{$last}{$hasha{$t}[1]-$i-$j}})
                                        {
                                            $tmp_multiDup=1; $str=1;
                                        }
                                        elsif($g-$j>0 and ($hasha{$t}[1]-$i-$j)>0 and exists $geneb{$k}{$g+$j} and exists $genea{$last}{$hasha{$t}[1]-$j} and exists $genea{$last}{$hasha{$t}[1]-$i-$j} and exists $pairs{$geneb{$k}{$g+$j}}{$genea{$last}{$hasha{$t}[1]-$j}} and exists $pairs{$geneb{$k}{$g+$j}}{$genea{$last}{$hasha{$t}[1]-$i-$j}})
                                        {
                                            $tmp_multiDup=1; $str=2;
                                        }
                                        else{last;}
                                    }else{
                                        if($str==1){
                                            unless($g-$j>0 and ($hasha{$t}[1]-$i-$j)>0 and exists $geneb{$k}{$g-$j} and exists $genea{$last}{$hasha{$t}[1]-$j} and exists $genea{$last}{$hasha{$t}[1]-$i-$j} and exists $pairs{$geneb{$k}{$g-$j}}{$genea{$last}{$hasha{$t}[1]-$j}} and exists $pairs{$geneb{$k}{$g-$j}}{$genea{$last}{$hasha{$t}[1]-$i-$j}})
                                            {$tmp_multiDup=0;}
                                        }elsif($str==2){
                                            unless($g-$j>0 and ($hasha{$t}[1]-$i-$j)>0 and exists $geneb{$k}{$g+$j} and exists $genea{$last}{$hasha{$t}[1]-$j} and exists $genea{$last}{$hasha{$t}[1]-$i-$j} and exists $pairs{$geneb{$k}{$g+$j}}{$genea{$last}{$hasha{$t}[1]-$j}} and exists $pairs{$geneb{$k}{$g+$j}}{$genea{$last}{$hasha{$t}[1]-$i-$j}})
                                            {$tmp_multiDup=0;}
                                        }
                                    }
                                }
                                if($tmp_multiDup==1){$tagDb{$geneb{$k}{$g}}=$i;}
                            }else{
                                $tagDb{$geneb{$k}{$g}}=1;
                            }
                        }elsif($last eq $hasha{$t}[0] and exists $genea{$last}{$hasha{$t}[1]+$i} and exists $pairs{$geneb{$k}{$g}}{$genea{$last}{$hasha{$t}[1]+$i}})
                        {
                            if($i>1){
                                my $tmp_multiDup=0;
                                my $str;
                                for(my $j=1; $j<=$i-1; $j++){
                                    if($j==1){
                                        if(exists $geneb{$k}{$g+$j} and exists $genea{$last}{$hasha{$t}[1]+$j} and exists $genea{$last}{$hasha{$t}[1]+$i+$j} and exists $pairs{$geneb{$k}{$g+$j}}{$genea{$last}{$hasha{$t}[1]+$j}} and exists $pairs{$geneb{$k}{$g+$j}}{$genea{$last}{$hasha{$t}[1]+$i+$j}})
                                        {
                                            $tmp_multiDup=1;$str=1;
                                        }
                                        elsif($g-$j>0 and exists $geneb{$k}{$g-$j} and exists $genea{$last}{$hasha{$t}[1]+$j} and exists $genea{$last}{$hasha{$t}[1]+$i+$j} and exists $pairs{$geneb{$k}{$g-$j}}{$genea{$last}{$hasha{$t}[1]+$j}} and exists $pairs{$geneb{$k}{$g-$j}}{$genea{$last}{$hasha{$t}[1]+$i+$j}})
                                        {
                                            $tmp_multiDup=1;$str=2;
                                        }else{last;}
                                    }else{
                                        if($str==1){
                                            unless(exists $geneb{$k}{$g+$j} and exists $genea{$last}{$hasha{$t}[1]+$j} and exists $genea{$last}{$hasha{$t}[1]+$i+$j} and exists $pairs{$geneb{$k}{$g+$j}}{$genea{$last}{$hasha{$t}[1]+$j}} and exists $pairs{$geneb{$k}{$g+$j}}{$genea{$last}{$hasha{$t}[1]+$i+$j}})
                                            {$tmp_multiDup=0;}
                                        }elsif($str==2){
                                            unless($g-$j>0 and exists $geneb{$k}{$g-$j} and exists $genea{$last}{$hasha{$t}[1]+$j} and exists $genea{$last}{$hasha{$t}[1]+$i+$j} and exists $pairs{$geneb{$k}{$g-$j}}{$genea{$last}{$hasha{$t}[1]+$j}} and exists $pairs{$geneb{$k}{$g-$j}}{$genea{$last}{$hasha{$t}[1]+$i+$j}})
                                            {$tmp_multiDup=0;}
                                        }
                                    }
                                }
                                if($tmp_multiDup==1){$tagDb{$geneb{$k}{$g}}=$i;}
                            }else{
                                $tagDb{$geneb{$k}{$g}}=1;
                            }
                        }else{
                            $tagDb{$geneb{$k}{$g}}=0,unless($tagDb{$geneb{$k}{$g}}>=1);
                            $last=$hasha{$t}[0];
                        }
                    }
###cycle for tendam length
                }
            }
        }
    }
}

##synteny score: 100 ~ 40
foreach my $ge (sort keys %hasha){
    foreach my $ge2 (keys %{$pairs{$ge}}){
        my @a = split(/\t/,$pairs{$ge}{$ge2});

        if($tagSa{$ge}==1 and $tagSb{$ge2}==1){
            if($a[0] eq $a[1]){
                print "$ge\t$tagSa{$ge}\t$tagSb{$ge2}\t$a[3]\t$a[2]\t+\t$ge2\t$tagDa{$ge}\t$tagDb{$ge2}\tSingleCopy\t$a[2]\t100\n";
                print "$ge2\t$tagSb{$ge2}\t$tagSa{$ge}\t$a[3]\t$a[2]\t+\t$ge\t$tagDb{$ge2}\t$tagDa{$ge}\tSingleCopy\t$a[2]\t100\n";
            }else{
                print "$ge\t$tagSa{$ge}\t$tagSb{$ge2}\t$a[3]\t$a[2]\t-\t$ge2\t$tagDa{$ge}\t$tagDb{$ge2}\tSingleCopy\t$a[2]\t100\n";
                print "$ge2\t$tagSb{$ge2}\t$tagSa{$ge}\t$a[3]\t$a[2]\t-\t$ge\t$tagDb{$ge2}\t$tagDa{$ge}\tSingleCopy\t$a[2]\t100\n";
            }
        }
        elsif($tagDa{$ge}>=1 or $tagDb{$ge2}>=1){
            if($a[0] eq $a[1]){
                print "$ge\t$tagSa{$ge}\t$tagSb{$ge2}\t$a[3]\t$a[2]\t+\t$ge2\t$tagDa{$ge}\t$tagDb{$ge2}\tTandem\t$a[2]\t40\n";
                print "$ge2\t$tagSb{$ge2}\t$tagSa{$ge}\t$a[3]\t$a[2]\t+\t$ge\t$tagDb{$ge2}\t$tagDa{$ge}\tTandem\t$a[2]\t40\n";
            }else{
                print "$ge\t$tagSa{$ge}\t$tagSb{$ge2}\t$a[3]\t$a[2]\t-\t$ge2\t$tagDa{$ge}\t$tagDb{$ge2}\tTandem\t$a[2]\t40\n";
                print "$ge2\t$tagSb{$ge2}\t$tagSa{$ge}\t$a[3]\t$a[2]\t-\t$ge\t$tagDb{$ge2}\t$tagDa{$ge}\tTandem\t$a[2]\t40\n";
            }
        }elsif($a[3] eq "Orphan"){
            if($a[0] eq $a[1]){
                print "$ge\t$tagSa{$ge}\t$tagSb{$ge2}\t$a[3]\t$a[2]\t+\t$ge2\t$tagDa{$ge}\t$tagDb{$ge2}\tOrphan\t$a[2]\t50\n";
                print "$ge2\t$tagSb{$ge2}\t$tagSa{$ge}\t$a[3]\t$a[2]\t+\t$ge\t$tagDb{$ge2}\t$tagDa{$ge}\tOrphan\t$a[2]\t50\n";
            }else{
                print "$ge\t$tagSa{$ge}\t$tagSb{$ge2}\t$a[3]\t$a[2]\t-\t$ge2\t$tagDa{$ge}\t$tagDb{$ge2}\tOrphan\t$a[2]\t50\n";
                print "$ge2\t$tagSb{$ge2}\t$tagSa{$ge}\t$a[3]\t$a[2]\t-\t$ge\t$tagDb{$ge2}\t$tagDa{$ge}\tOrphan\t$a[2]\t50\n";
            }
        }else{
            my $tmp_score=60;
            if($a[3] eq "Perfect"){$tmp_score=90;}
            elsif($a[3] eq "Good"){$tmp_score=80;}
            elsif($a[3] eq "OK"){$tmp_score=70;}
            elsif($a[3] eq "UpNone" or $a[3] eq "DownNone"){$tmp_score=60;}

            if($a[0] eq $a[1]){
                print "$ge\t$tagSa{$ge}\t$tagSb{$ge2}\t$a[3]\t$a[2]\t+\t$ge2\t$tagDa{$ge}\t$tagDb{$ge2}\tOthers\t$a[2]\t$tmp_score\n";
                print "$ge2\t$tagSb{$ge2}\t$tagSa{$ge}\t$a[3]\t$a[2]\t+\t$ge\t$tagDb{$ge2}\t$tagDa{$ge}\tOthers\t$a[2]\t$tmp_score\n";
            }else{
                print "$ge\t$tagSa{$ge}\t$tagSb{$ge2}\t$a[3]\t$a[2]\t-\t$ge2\t$tagDa{$ge}\t$tagDb{$ge2}\tOthers\t$a[2]\t$tmp_score\n";
                print "$ge2\t$tagSb{$ge2}\t$tagSa{$ge}\t$a[3]\t$a[2]\t-\t$ge\t$tagDb{$ge2}\t$tagDa{$ge}\tOthers\t$a[2]\t$tmp_score\n";
            }
        }
    }
}
