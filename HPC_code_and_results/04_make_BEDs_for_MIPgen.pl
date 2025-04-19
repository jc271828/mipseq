#!/usr/bin/perl -w

use strict;

if ($#ARGV !=4) {
    die "usage: $0 SNPs flank skipEnds refLength outSuffix\n";
}

my %data=();
open (IN1, "<$ARGV[0]") || die $!;
while (<IN1>) {
    chomp;
    my @tmp=split(/\t/);
    push (@{$data{$tmp[4]}}, [@tmp]);
}
warn scalar keys %data, " strains read\n";

my %ref=();
open (IN2, "<$ARGV[3]") || die $!;
while (<IN2>) {
    chomp;
    my @tmp=split(/\t/);
    $ref{$tmp[0]}=$tmp[1];
}
warn scalar keys %ref, " chromosomes read\n";


foreach my $strain (sort {$a cmp $b} keys %data) {
    my $fname=$strain.$ARGV[4];
    open(OUT, ">$fname") || die $!;
    foreach (@{$data{$strain}}) {
	my @line=@{$_};
	if ($line[1]-$ARGV[2] < 1) {
	    next;
	}
	if ($line[1]+$ARGV[2] > $ref{$line[0]}) {
	    next;
	}
	my $start=$line[1]-$ARGV[1];
	my $stop=$line[1]+$ARGV[1];
	my $id=$line[4]."_$line[0]:$line[1]:$line[2]->$line[3]";
	print OUT "$line[0]\t$start\t$stop\t$id\n";
    }
    close OUT;
}

