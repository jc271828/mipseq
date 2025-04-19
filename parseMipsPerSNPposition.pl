#!/usr/bin/perl -w

use strict;

if ($#ARGV != 1) {
    die "usage: $0 mips distanceFromLigProbeEnd\n"; # distance includes lig probe length, but not UMI, which is by default at the extension probe side
}

open (IN, "<$ARGV[0]") || die "cannot open $ARGV[0]: $!\n";

my $total=0;
my $count=0;

while (<IN>) {
    chomp;
    my @tmp=split(/\t/);
    if ($tmp[0]=~/mip_key/) {
	print "$_\tmin_read_length\n";
	next;
    }
    $total++;
    my $lig_probe_length=length($tmp[10]);
    my $strand=$tmp[17];
    my $snp_pos=$tmp[19]=~/\S+:(\d+):\S+/ ? $1 : die "cannot parse SNP position $tmp[19]\n";
#    print "$lig_probe_length\t$strand\t$snp_pos\n";
    my $dist_to_snp;
    if ($strand eq "+") {
	$dist_to_snp=$tmp[8]-$snp_pos+1; # lig_probe_stop
    }
    else {
	$dist_to_snp=$snp_pos-$tmp[7]+1; # lig_probe_start
    }
    if ($dist_to_snp <= $ARGV[1] && $dist_to_snp >= $lig_probe_length) { # scan sequence may not cover the SNP position
	print join("\t", @tmp), "\t$dist_to_snp\n";
	$count++;
    }    
}

warn "$total MIPs read\n";
warn "$count MIPs have SNPs within $ARGV[1] bp from start of read 1\n";
