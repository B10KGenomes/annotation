#!usr/bin/perl -w
use strict;

open(IN,"$ARGV[0]") or die $!;
$/=">";<IN>;$/="\n";
while (<IN>) {
	chomp;
	$_ =~ /(\S+)\s/;
	my $id=$1;
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$/="\n";
	$seq =~ s/\n//g;
	my $t_len = length ($seq);
	$seq =~ s/[A-Z]//g;
	my $m_len = length ($seq);
	my $rate = $m_len/$t_len;
	print "$id\t$rate\n";
}
close IN;
