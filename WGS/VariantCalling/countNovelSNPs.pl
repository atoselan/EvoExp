#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw( min max );
use List::MoreUtils qw(uniq);

#--- Count novel SNPs
# Count number of novel SNPs in each sample
# where T0 (first 3 samples) is homozygous to
# the reference but is either heterozygous or
# homozygous to the alternate base in the
# evolved sample.

# Use the following command to filter bcf file first:
# bcftools view -V indels -M2
# This removes indels and retains only bi-allelic SNPs


my ($vcf) = @ARGV;


### Set filtering cutoffs
my $minQual       = 40;         # minimum variant quality
my $minDepth      = 10;         # minimum read depth
my $maf           = 0.05;       # Minimum minor allele frequency (MAF) to consider a het snp real
my $minHomSupport = 0.95;       # Required read support for hom snp
my $maxOccur      = 5;          # If var occurs multiple times, more likely to be a rare allele missed in T0
my %novelCount;


open (IN, "<$vcf") or die "Can't open $vcf: $!\n";
while (<IN>)
{
  my $line = $_;
  chomp $line;
  my $hetCount = 0;
  my $homCount = 0;
  # If not a comment line
  if ($line !~ m/^#/)
  {
    $line =~ s/\s+/\t/g;                # Replace multi space with tab
    my @fields = split(/\t/,$line);     # Split line on tab symbol
    # Check variant has passed filters and qual score is greater than threshold
    if (($fields[6] eq 'PASS') && ($fields[5] >= $minQual))
    {
      ### Get control sample info. For some reason there are two different formats for this... Check for both.
      # Get consensus genotype
      my $consensusGenotype = &getConsensusGenotype(\@fields);
      if ($consensusGenotype ne 'FAIL')
      {
        #print "$consensusGenotype";
        # Loop through individual samples starting with sample after control
        # count how many times SNP occurs
        for (my $i=12; $i <= $#fields; $i++)
        {
          # Get format fields for sample
          my ($gt,$pl,$dp,$ad)=0;
          if ($fields[8] eq 'GT:PL:DP:AD')
          {
            ($gt,$pl,$dp,$ad) = split(/:/,$fields[$i]);
          }
          elsif ($fields[8] eq 'GT:DP:AD:PL')
          {
            ($gt,$dp,$ad,$pl) = split(/:/,$fields[$i]);
          }
          #print "\t$gt";
          if ($gt !~ m/\./) #&& $ad !~ m/\./)
          {
            my @alleleDepths = split(/,/,$ad);
            if ($gt eq '0/0')
            {
              pop(@alleleDepths);
            }
            my $lowestAlleleDepth = min(@alleleDepths);
            my $highestAlleleDepth = max(@alleleDepths);
            # If control is homozygous:
            if (($consensusGenotype eq '1/1' || $consensusGenotype eq '0/0') && ($dp >= $minDepth))
            {
              # If evolved line is heterozygous:
              if (($gt eq '0/1' || $gt eq '0/2' || $gt eq '1/2') &&
                  $lowestAlleleDepth && (($lowestAlleleDepth / $dp) >= $maf))
              {
                $hetCount++;
              }
              elsif ((($consensusGenotype eq '0/0' && $gt eq '1/1' ) || ($consensusGenotype eq '1/1' && $gt eq '0/0')) &&
                     $highestAlleleDepth && (($highestAlleleDepth / $dp) >= $minHomSupport))
              {
                $homCount++;
              }
            }
          }
        }
        #print "$homCount\t$hetCount\n";
        # Loop thorugh again if number of occurrences of SNP is less than maxOccur threshold
        my $snpFound = 0;
        if (($consensusGenotype eq '0/0' || $consensusGenotype eq '1/1') && ($hetCount + $homCount) <= $maxOccur)
        {
          #print "1\n";
          #my $homRflag = 0;
          for (my $i=12; $i <= $#fields; $i++)
          {
            my ($gt,$pl,$dp,$ad) = '';
            if ($fields[8] eq 'GT:PL:DP:AD')
            {
              ($gt,$pl,$dp,$ad) = split(/:/,$fields[$i]);
              #print "$gt\n";
            }
            elsif ($fields[8] eq 'GT:DP:AD:PL')
            {
              ($gt,$dp,$ad,$pl) = split(/:/,$fields[$i]);
              #print "$gt\n";
            }
            if ($gt !~ m/\./) #&& $ad !~ m/\./)
            {
              #print "2\n"  ;
              my @alleleDepths = split(/,/,$ad);
              if ($gt eq '0/0')
              {
                pop(@alleleDepths);
              }
              my $lowestAlleleDepth = min(@alleleDepths);
              my $highestAlleleDepth = max(@alleleDepths);
              if ($gt ne './.')
              {
                if ((($gt eq '0/1' || $gt eq '0/2' || $gt eq '1/2') &&
                    $lowestAlleleDepth &&
                    ($dp >= $minDepth) &&
                    (($lowestAlleleDepth / $dp) >= $maf)) ||
                    (($consensusGenotype eq '0/0' && $gt eq '1/1') && $highestAlleleDepth && (($highestAlleleDepth / $dp) >= $minHomSupport)) ||
                    (($consensusGenotype eq '1/1' && $gt eq '0/0') && $highestAlleleDepth && (($highestAlleleDepth / $dp) >= $minHomSupport)))
                {
                  $novelCount{$i}++;
                  $snpFound = 1;
                  #print "$i,";
                }
              }
            }
          }
          # EDIT: changed this line - removed tab at start of line
          if ($snpFound){print "$line\n";}
        }
      }
    }
  }
  # Else, print comment line as is
  else
  {
    print "$line\n";
  }
} # while
close IN or die "Can't close $vcf: $!\n";



# Print total novel SNPs per sample
#foreach my $key (sort keys %novelCount)
#{
#  print "$key\t$novelCount{$key}\n";
#}



# Check if control sample genotype is consistent
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
    if ($gt eq '0/0' || $gt eq '1/1' && $dp >= $minDepth)
    {
      # Get highest and lowest AF
      my @alleleDepths = split(/,/,$ad);
      if ($gt eq '0/0')
      {
        pop(@alleleDepths);
      }
      my $lowestAlleleDepth = min(@alleleDepths);
      my $highestAlleleDepth = max(@alleleDepths);
      # Calculate support for homozygous position
      my $homSupport = 0;
      if ($highestAlleleDepth && $dp)
      {
        $homSupport = $highestAlleleDepth / $dp;
      }
      # If min read support and homo support met push genotype onto gt array
      if ($dp >= $minDepth && $homSupport >= $minHomSupport)
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
  if (scalar(uniq(@genotypes)) == 1 && ($genotypes[0] eq '0/0' || $genotypes[0] eq '1/1'))
  {
    return $genotypes[0];
  }
  else
  {
    return "FAIL";
  }
}






=pod

Change script to count number of novel SNPs

Have tweaked the SNP calling script to include DP and AD in FORMAT string (GT:PL:DP:AD)

Have added two extra replicates of T0. Need to compare to all 3 rather than just one.

  - Get consensus of the 3 T0 samples.
    - If genotype is the same for all 3
    - If depth is ok for all 3
    - If MAF is ok for all 3


# Use a 3 element array for gt
@gt


for (1st to 3rd sample)
  get gt,dp,ad,pl
  get high and lowest AF
  hom_support = highest_AF / depth
  if (dp>= min_depth && hom_support >= min_hom)
    push gt onto array



if ((length uniq(@gt)) == 1 and (first element either 0/0 or 1/1))
  then consistent homozygous call
  store first element as consensus for this row

  # Check all other samples against this

Also, can I set the looping through the other samples as a subroutine?
I do this twice, so makes sense.

=cut