#!/usr/bin/perl
use strict;
use warnings;

#--- formatWeirFSTforCircos.pl
# Format a windowed FST results file for use with Circos. 
# Basically just involves using a new chromosome number and
# taking a subset of fields present in .fst file.


my ($karyoFile, $fstFile) = @ARGV;

my $minFST = 0;
my %cHash = &prepChrHash($karyoFile);


open (IN, "<$fstFile") or die "Can't open $fstFile: $!\n";
while(<IN>)
{
  my $line = $_;
  chomp $line;

  if ($line !~ m/^CHROM/)
  {
    my ($chr,$start,$end,$nVar,$wFST,$mFST) = split(/\t/, $line);
    # If greater than zero
    if ($wFST > $minFST)
    {
      # Check if end of window > chromosome length
      my $circId = $cHash{$chr}{circos};
      print "$circId $start ";
      if ($end > $cHash{$chr}{length})
      {
        print "" . $cHash{$chr}{length} . " ";
      }
      else
      {
        print "$end ";
      }
      print "$wFST\n";
    }
  }
}
close IN or die;



# Parse karyotype file. Use this to replace ID from originak
# fst file to labels for Circos. 
sub prepChrHash
{
  my $file = shift;
  my %chrHash;
  open (IN, "<$file") or die "Can't open $file: $!\n";
  while (<IN>)
  {
    my $line = $_;
    chomp $line;
    my ($prefix, $foo, $circosId, $chrId, $start, $length, $color) = split(/\s/, $line);
    # Create ID compatible with original chromosome IDs for lookup
    # Fc data has scaffold prefix, Tp starts with chr
    #my $id = $prefix . "_" . $chrId;
    my $id = "scaffold_" . $chrId;
    $chrHash{$id}{circos} = $circosId;
    $chrHash{$id}{length} = $length;
  }
  close IN or die;
  return %chrHash;
}





=pod

Vcftools Fst format
CHROM	BIN_START	BIN_END	N_VARIANTS	WEIGHTED_FST	MEAN_FST
chr_1	1	10000	275	0.00387777	-0.00754446
chr_1	1001	11000	279	0.0129056	0.00119948
chr_1	2001	12000	256	0.0183706	0.00567678
chr_1	3001	13000	239	0.0202634	0.00664598
chr_1	4001	14000	211	0.020566	0.00414267
chr_1	5001	15000	230	0.0197341	0.00286338
chr_1	6001	16000	217	0.0203849	0.00241217


Output for Circos
chr - tp1 1 0 3042585 chr1
chr - tp2 2 0 2707195 chr2
chr - tp3 3 0 2440052 chr3
chr - tp4 4 0 2402323 chr4
chr - tp5 5 0 2305972 chr5
chr - tp6 6 0 2071480 chr6
chr - tp7 7 0 1992434 chr7
chr - tp8 8 0 1267198 chr8
=cut
