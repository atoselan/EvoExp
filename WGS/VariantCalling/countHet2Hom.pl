#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw( min max );
use List::MoreUtils qw(uniq);

#--- Count het to hom
# Slight variation on previous script to count number of novel SNPs.
# This script counts the number of heterozygous SNPs in T0 that reach
# fixation in an evolved cell line.
#
# Use bcftools to format vcf file to feed into this script.
# E.g.
#
# bcftools view -M 2 -V indels -O v \
#   mergedFiltered_NoMismatchFilter_allBatches_SNPsOnly.vcf \
#   | perl countHet2Hom_allBatches.pl -


# Read in vcf file and optionally a specific chromosome
my ($vcf, $chr) = @ARGV;


### Set filtering cutoffs
my $minQual    = 40;    # minimum variant quality
my $minDepth   = 10;    # minimum read depth
my $minMaf     = 0.05;  # Minimum minor allele frequency (MAF) to consider a het snp real
my $homSupport = 0.95;  # Required read support for hom snp
my $maxOccur   = 35;    # If var occurs multiple times, more likely to be a rare allele missed in T0
my %fixedCount;


open (IN, "<$vcf") or die "Can't open $vcf: $!\n";
while (<IN>)
{
  my $line = $_;
  chomp $line;
  my $homCount = 0;     # Count how many samples SNP occurs in

  # If not a comment line
  if ($line !~ m/^#/)
  {
    $line =~ s/\s+/\t/g;                # Replace multi space with tab
    my @fields = split(/\t/,$line);     # Split line on tab symbol
    # Check variant has passed filters and qual score is greater than threshold
    if (($chr && ($fields[0] eq $chr)) || !$chr)
    {
      if ($fields[6] eq 'PASS' && $fields[5] >= $minQual)
      {
        ### Get control sample info. For some reason there are two different formats for this... Check for both.
        # Get consensus genotype
        my $consensusGenotype = &getConsensusGenotype(\@fields);
        if ($consensusGenotype ne 'FAIL')
        {
          # Loop through individual samples starting with sample after control
          # if SNP is homozygous in current sample, increment counter for that
          # sample. Use loop counter as index for hash.
          for (my $i=12; $i <= $#fields; $i++)
          {
            # Get format fields for sample
            my ($gt,$pl,$dp,$ad)='';
            if ($fields[8] eq 'GT:PL:DP:AD')
            {
              ($gt,$pl,$dp,$ad) = split(/:/,$fields[$i]);
            }
            elsif ($fields[8] eq 'GT:DP:AD:PL')
            {
              ($gt,$dp,$ad,$pl) = split(/:/,$fields[$i]);
            }

            #---
            #--- HERE
            #---
            # Need to change this!!!!
            # Error message about using . as numeric. Caused by missing genotype.
            if ($gt !~ m/\./ && ($dp >= $minDepth))
            {
              my @alleleDepths = split(/,/,$ad);
              if (!$alleleDepths[0]){$alleleDepths[0]=0;}
              if (!$alleleDepths[1]){$alleleDepths[1]=0;}
              my ($lowestAlleleDepth, $highestAlleleDepth) = '';
              if ($alleleDepths[0] eq "\."){$alleleDepths[0] = 0;}
              if ($alleleDepths[1] eq "\."){$alleleDepths[1] = 0;}
              $lowestAlleleDepth = min(@alleleDepths);
              $highestAlleleDepth = max(@alleleDepths);
              # If fixed meets min depth and hom support count number of occurrences
              if (($gt eq '1/1' || $gt eq '0/0') && (($highestAlleleDepth / $dp) >= $homSupport)) #&& ($dp >= $minDepth))
              {
                $homCount++;
              }
            }
          }

          # Hash to store the names of samples that have a particular fixed heterozygous SNP
          my %sampleHash;
          my $fixedFlag = 0;
          # Second iteration, if most are fixed hom ref we've probably just missed a rare allele
          # if most are fixed in the evolved we've probably missed it in T0
          if ($homCount && ($homCount <= $maxOccur))
          {
            #my $fixedFlag = 0;
            for (my $i=12; $i <= $#fields; $i++)
            {
              my ($gt,$pl,$dp,$ad) = '';
              if ($fields[8] eq 'GT:PL:DP:AD')
              {
                ($gt,$pl,$dp,$ad) = split(/:/,$fields[$i]);
              }
              elsif ($fields[8] eq 'GT:DP:AD:PL')
              {
                ($gt,$dp,$ad,$pl) = split(/:/,$fields[$i]);
              }
              if ($gt !~ m/\./ && ($dp >= $minDepth)) #&& $ad !~ m/\./)
              {
                my @alleleDepths = split(/,/,$ad);
                if (!$alleleDepths[0]){$alleleDepths[0]=0;}
                if (!$alleleDepths[1]){$alleleDepths[1]=0;}
                my ($lowestAlleleDepth, $highestAlleleDepth) = '';
                if ($alleleDepths[0] eq "\."){$alleleDepths[0] = 0;}
                if ($alleleDepths[1] eq "\."){$alleleDepths[1] = 0;}
                $lowestAlleleDepth = min(@alleleDepths);
                $highestAlleleDepth = max(@alleleDepths);
                if (($gt eq '1/1' || $gt eq '0/0') && (($highestAlleleDepth / $dp) >= $homSupport)) #&& ($dp >= $minDepth))
                {
                  $fixedCount{$i}++;
                  $fixedFlag = 1;
                  $sampleHash{$i} = 1;
                }
              }
            }
          }
          if ($fixedFlag)
          {
            print "$line\n";
            #print "$fields[0]\t$fields[1]";
            #for (my $i=10; $i <= $#fields; $i++)
            #{
            #  (exists $sampleHash{$i})? print "\t1": print "\t0";
            #}
            #print "\n";
          }
        }
      }
    }
  }
  # Else if comment line print
  else
  {
    print "$line\n";
  }
} # while
close IN or die "Can't close $vcf: $!\n";


#foreach my $key (sort keys %fixedCount)
#{
#  print "$key\t$fixedCount{$key}\n";
#}


# Is there a consensus between the genotype of all three T0
# samples? Do they meet the specified depth criteria?
sub getConsensusGenotype
{
  my $refFields = $_[0];
  my @fields = @$refFields;
  my @genotypes;
  my $startIndx = 9;
  my $endIndx = 11;
  # Loop through elements 9 to 11 (these are the three control samples)
  for (my $i=$startIndx; $i<=$endIndx; $i++)
  {
    # Check format of FORMAT field, there are two (that I can see) different ones
    my ($gt,$pl,$dp,$ad) = 0;
    if ($fields[8] eq 'GT:PL:DP:AD')
    {
      ($gt,$pl,$dp,$ad) = split(/:/, $fields[$i]);
    }
    elsif ($fields[8] eq 'GT:DP:AD:PL')
    {
      ($gt,$dp,$ad,$pl) = split(/:/, $fields[$i]);
    }
    # If heterozygous
    if ($gt !~ m/\./ && ($gt eq '0/1' || $gt eq '0/2' || $gt eq '1/2') && $dp >= $minDepth)
    {
      # Get highest and lowest AF
      my @alleleDepths = split(/,/,$ad);
      if (!$alleleDepths[0]){$alleleDepths[0]=0;}
      if (!$alleleDepths[1]){$alleleDepths[1]=0;}
      my ($lowestAlleleDepth, $highestAlleleDepth) = '';
      if ($alleleDepths[0] =~ m/\./){$alleleDepths[0] = 0;}
      if ($alleleDepths[1] =~ m/\./){$alleleDepths[1] = 0;}
      #print "@alleleDepths\n";
      $lowestAlleleDepth  = min(@alleleDepths);
      $highestAlleleDepth = max(@alleleDepths);
      # Calculate MAF
      my $maf = $lowestAlleleDepth / $dp;
      # If min read support and MAF thresholds met push genotype onto gt array
      if ($maf >= $minMaf)
      {
        push(@genotypes,$gt);
      }
      else
      {
        push(@genotypes,"FAIL");
      }
    }
    else
    {
      push(@genotypes,"FAIL");
    }
  }
  if (scalar(@genotypes ==3) && scalar(uniq(@genotypes)) == 1 && ($genotypes[0] eq '0/1' || $genotypes[0] eq '0/2' || $genotypes[0] eq '1/2'))
  {
    return $genotypes[0];
  }
  else
  {
    #else return a don't care
    return "FAIL";
  }
}