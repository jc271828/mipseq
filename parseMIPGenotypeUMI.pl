#!/usr/bin/perl -w

use strict;

if ($#ARGV!=2) {
    die "usage: $0 MIPs reads UMIlength\n";
}

open (IN1, "<$ARGV[0]") || die $!;
if ($ARGV[1]=~/\.gz/) {
    open (IN2, '-|', "zcat $ARGV[1]") || die $!;
}
else {
  open (IN2, "<$ARGV[1]") || die $!;
}

my $sample_id=$ARGV[1]=~/^(\S+)\.fastq\.gz/ ? $1 : die " cannot parse sample ID: $ARGV[1]\n";
warn "sample ID: $sample_id\n";

my %mips=();
while(<IN1>) {
  chomp;
  my @tmp=split(/\t/);
  my  ($strain, $chr, $pos, $ref, $alt)=$tmp[19]=~/^(\S)+_(\S+):(\d+):(\S)->(\S)/ ? ($1, $2, $3, $4, $5) : die "cannot parse SNP $tmp[19]\n";
  %{$mips{$tmp[19]}} = ("strain" => $strain,
			"chr" => $chr,
			"pos" => $pos,
			"ref" => $ref,
			"alt" => $alt
		       );
  $mips{$tmp[19]}{"strand"}=$tmp[17];
  my $revcompligprobe=reverse($tmp[10]);
  $revcompligprobe=~tr/ACGT/TGCA/;
  $mips{$tmp[19]}{"revcompligprobe"}=$revcompligprobe; 
  $mips{$tmp[19]}{"data"}=\@tmp;

  my $snpPos;
  if ($mips{$tmp[19]}{strand} eq '+') {
    $snpPos=$mips{$tmp[19]}{pos}-$mips{$tmp[19]}{data}[11]+0;  #0-based; mip_scan_start_position
  }
  else {
    $snpPos=$mips{$tmp[19]}{data}[12]-$mips{$tmp[19]}{pos}+0;  #0-based; mip_scan_stop_position
  }
  my @target=split("|", $mips{$tmp[19]}{data}[13]);
  my $seq='';
  my $revComp='';
  if ($snpPos+1 <= length($mips{$tmp[19]}{data}[13])) {      
      $seq=substr($mips{$tmp[19]}{data}[13], $snpPos+1);  #not including SNP
      $revComp=reverse $seq;
      $revComp=~tr/AGCT/TCGA/;
  }
  else {
      warn "$tmp[0]\t$tmp[19]\t", $snpPos+1, "\t$mips{$tmp[19]}{data}[13]\t", length($mips{$tmp[19]}{data}[13]), "\n";
  }
  my $revCompReadUpToSnp=$mips{$tmp[19]}{revcompligprobe}.$revComp;
  my $snpPosInRead=length($revCompReadUpToSnp);
  $mips{$tmp[19]}{"revCompReadUpToSnp"}=$revCompReadUpToSnp;
  $mips{$tmp[19]}{"snpPosInRead"}=$snpPosInRead;

  $mips{$tmp[19]}{refCount}=0;
  $mips{$tmp[19]}{altCount}=0;
  %{$mips{$tmp[19]}{refUmi}}=();
  %{$mips{$tmp[19]}{altfUmi}}=();
=head
  print "$tmp[19]\n";
  print "$mips{$tmp[19]}{strand}\n";
  print "$mips{$tmp[19]}{data}[11]\t$mips{$tmp[19]}{data}[12]\n";
  print "$mips{$tmp[19]}{data}[13]\n";
  print join("-", @target), "\n";
  print "$snpPos\t$target[$snpPos]\n";
  print "$seq\n";
  print "$revComp\n";
  print "$mips{$tmp[19]}{revcompligprobe}\n";
  print "$revCompReadUpToSnp\n";
=cut
}
warn scalar keys %mips, " MIPs read\n";

		       


my $i=0;
my $readId='';
my $readCount=0;
my $withMips=0;
my $noMips=0;
my $seqMatch=0;

my %umiCounts=();
    
while (<IN2>) {
  chomp;
  $i++;
  if ($i==1) {
    $readId=$_=~/\@(\S+)/ ? $1 : die "cannot parse read ID: $_";
    $readCount++;
  }
  elsif ($i==2) {
    my $umi=substr($_, 0, $ARGV[2]);
    my $seq=substr($_, $ARGV[2]);
    my $flag=0;
    $umiCounts{total}{$umi.$seq}++;
    foreach my $snp (sort {$a cmp $b} keys %mips) {
      if ($seq=~/^$mips{$snp}{"revcompligprobe"}/) {
	$withMips++;
	$flag=1;
	$umiCounts{withMips}{$umi.$mips{$snp}{"revcompligprobe"}}++;
	if ($seq=~/^$mips{$snp}{"revCompReadUpToSnp"}/) {  # match lig_probe and up to SNP in scan sequence
	    $seqMatch++;
	    $umiCounts{seqMatch}{$umi.$mips{$snp}{"revCompReadUpToSnp"}}++;
	  my $base=substr($seq, $mips{$snp}{snpPosInRead}, 1);
	  #	  warn "$seq\n$mips{$snp}{revcompligprobe}\n$mips{$snp}{revCompReadUpToSnp}\n$base\n", substr($seq, $mips{$snp}{snpPosInRead}, 5), "\n"; # correct
	  if ($mips{$snp}{strand} eq '+') {
	    $base=~tr/AGCT/TCGA/;
	  }
	  if ($base eq $mips{$snp}{ref}) {
	    $mips{$snp}{refCount}++;
	    $mips{$snp}{refUmi}{$umi}++;
	  }
	  elsif ($base eq $mips{$snp}{alt}) {
	    $mips{$snp}{altCount}++;
	    $mips{$snp}{altUmi}{$umi}++;
	  }
	  else {	    
	    $mips{$snp}{other}{$base}{count}++;
	    $mips{$snp}{other}{$base}{umi}{$umi}++;
#	    warn "$snp\t$mips{$snp}{strand}\n$seq\n$mips{$snp}{revcompligprobe}\n$mips{$snp}{revCompReadUpToSnp}\n$base\n", substr($seq, $mips{$snp}{snpPosInRead}, 5), "\n";
	  }
	}
	else {
#	  warn "$seq\n$mips{$snp}{revcompligprobe}\n$mips{$snp}{revCompReadUpToSnp}\n"; #seemd correct - at least one mismatch
	}
	last;
      }
    }
    if (!$flag) {
	$noMips++;
	$umiCounts{noMips}{$umi.$seq}++;
    }
  }
  elsif ($i==4) {
    $i=0;
  }
}

my $umiCountsTotal=scalar keys %{$umiCounts{total}};
my $umiCountsWithMips=scalar keys %{$umiCounts{withMips}};
my $umiCountsSeqMatch=scalar keys %{$umiCounts{seqMatch}};
my $umiCountsNoMips=scalar keys %{$umiCounts{noMips}};

warn "$readCount reads read ($umiCountsTotal unique UMIs)\n";
warn "$withMips reads with MIPs identified ($umiCountsWithMips unique UMIs)\n";
warn "$seqMatch reads match sequence to SNP ($umiCountsSeqMatch unique UMIs)\n";
warn "$noMips reads had no corresponding MIPs ($umiCountsNoMips unique UMIs)\n";

print "# sample ID:\t$sample_id\n";
print "# reads in library:\t$readCount\t($umiCountsTotal unique UMIs)\n";
print "# reads with lig_probe sequence:\t$withMips\t($umiCountsWithMips unique UMIs)\n";
print "# reads with lig_probe and scan sequence (up to SNP):\t$seqMatch\t($umiCountsSeqMatch unique UMIs)\n";
print "# reads with no match:\t$noMips\t($umiCountsNoMips unique UMIs)\n";
print "mip_name\tref\talt\tref_read_count\talt_read_count\tother1_read_count\tother2_read_count\tref_read_freq\talt_read_freq\tother1_read_freq\tother2_read_freq\tref_umi_count\talt_umi_count\tother1_umi_count\tother2_umi_count\tref_umi_freq\talt_umi_freq\tother1_umi_freq\tother2_umi_freq\n";
foreach my $snp (sort {$a cmp $b} keys %mips) {
  my $refUmiCount=scalar keys %{$mips{$snp}{refUmi}};
  my $altUmiCount=scalar keys %{$mips{$snp}{altUmi}};
  print "$snp\t$mips{$snp}{ref}\t$mips{$snp}{alt}\t$mips{$snp}{refCount}\t$mips{$snp}{altCount}";
  #  foreach my $base (sort {$a cmp $b} keys %{$mips{$snp}{other}}) {
  foreach my $base ("A", "C", "G", "T") {
    if ($base eq $mips{$snp}{ref} || $base eq $mips{$snp}{alt}) {
      next;
    }
    else {
      print "\t$base:";
      if ($mips{$snp}{other}{$base}{count}) {
	print "$mips{$snp}{other}{$base}{count}";
      }
      else {
	print "0";
      }
    }
  }
  my ($ref_freq, $alt_freq, $other_freq)=(0,0,0);
  $ref_freq=sprintf("%.6f", $mips{$snp}{refCount}/$seqMatch);
  $alt_freq=sprintf("%.6f", $mips{$snp}{altCount}/$seqMatch);
  print "\t$ref_freq\t$alt_freq";
  foreach my $base ("A", "C", "G", "T") {
      if ($base eq $mips{$snp}{ref} || $base eq $mips{$snp}{alt}) {
	  next;
      }
      else {
	  print "\t$base:";
	  if ($mips{$snp}{other}{$base}{count}) {
	      $other_freq=sprintf("%.6f", $mips{$snp}{other}{$base}{count}/$seqMatch);
	      print "$other_freq";
	  }
	  else {
	      print "0";
	  }
      }
  }
#  print "\n";
  print "\t$refUmiCount\t$altUmiCount";
#  foreach my $base (sort {$a cmp $b} keys %{$mips{$snp}{other}}) {
  foreach my $base ("A", "C", "G", "T") {
    if ($base eq $mips{$snp}{ref} || $base eq $mips{$snp}{alt}) {
      next;
    }
    else {
      print "\t$base:";
      if ($mips{$snp}{other}{$base}{umi}) {
	my $umiCount=scalar keys %{$mips{$snp}{other}{$base}{umi}};
	print "$umiCount";
      }
      else {
	print "0";
      }
    }
  }
  $ref_freq=sprintf("%.6f", $refUmiCount/$umiCountsSeqMatch);
  $alt_freq=sprintf("%.6f", $altUmiCount/$umiCountsSeqMatch);
  print "\t$ref_freq\t$alt_freq";
  foreach my $base ("A", "C", "G", "T") {
      if ($base eq $mips{$snp}{ref} || $base eq $mips{$snp}{alt}) {
	  next;
      }
      else {
	  print "\t$base:";
	  if ($mips{$snp}{other}{$base}{umi}) {
	      my $umiCount=scalar keys %{$mips{$snp}{other}{$base}{umi}};
	      $other_freq=sprintf("%.6f", $umiCount/$umiCountsSeqMatch);
	      print "$other_freq"
	  }
	  else {
	      print "0";
	  }
      }
  }

  print "\n";
}

