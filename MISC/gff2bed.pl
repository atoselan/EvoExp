#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min max);

# Convert gff file to bed format
# One line per gene, need to determine the start and 
# stop positions of the gene overall, how many exons
# in total it is comprised of, the start and stop 
# postions of each exon.

# See the link below for description of Bed format:
# http://genome.ucsc.edu/FAQ/FAQformat


my $gff = shift;


# Hash to store values
my %prot_hash;
# Store current prot id
my $tmp_prot = '';
my $prot = '';

open (IN, "<$gff") or die "Can't open $gff: $!\n";
while (<IN>)
{
  my $line = $_;
  chomp $line;
  # If coding region
  if ($line =~ m/.*CDS.*/)
  {
    # Get fields in line
    my @fields = split(/\t/, $line);
    my $chr = $fields[0];
    # Remove underscore from chromosome name, it will crash methylkit
    $chr =~ s/\_//;
    my $start = $fields[3];
    my $stop = $fields[4];
    my $strand = $fields[6];
    my $details = $fields[8];

    if ($details =~ m/.*proteinId\s(\d+);/)
    {
      $prot = $1;
      #print "$prot\t$tmp_prot\n";
    }

    # If protId is different from tmp then we have a new value and
    # should print out the current value
    if (($prot && $tmp_prot) && ($prot ne $tmp_prot))
    {
      # print "$prot\t$tmp_prot\n";
      print "$prot_hash{$tmp_prot}{chrom}\t";		# Chromosome name
      # Determine gene start and stop positions
      my @starts = split(/,/, $prot_hash{$tmp_prot}{start});
      my @stops = split(/,/, $prot_hash{$tmp_prot}{stop});
      my $prot_start = min(@starts);
      my $prot_stop = max(@stops);
      print "$prot_start\t";				# Protein start
      print "$prot_stop\t";				# Protein end
      print "Thaps-$tmp_prot\t";			# Protein name
      print "0\t";					# Score - not needed
      print "$prot_hash{$tmp_prot}{strand}\t";		# Strand
      print "$prot_start\t";                            # Thick start
      print "$prot_stop\t";                             # Thick end
      print "0\t";					# itemRgb - not needed
      print "$prot_hash{$tmp_prot}{block_count}\t";	# Block count, num of exons
      print "$prot_hash{$tmp_prot}{size},\t";		# Block sizes
      # Get start values in relation to gene start
      foreach (@starts)
      {
        print "" . ($_ - $prot_start) . ",";
      }
      print "\n";
      #print "$prot_hash{$tmp_prot}{start},\n";		# Block starts
    }

    if (exists $prot_hash{$prot})
    {
      $prot_hash{$prot}{start} .= ",";
      $prot_hash{$prot}{stop} .= ",";
      $prot_hash{$prot}{size} .= ",";
      $prot_hash{$prot}{block_count}++;
    }
    # If new prot initialise some of the values
    else
    {
      $prot_hash{$prot}{chrom} = $chr; 
      $prot_hash{$prot}{block_count} = 1;
      $prot_hash{$prot}{strand} = $strand;
    }
    $prot_hash{$prot}{start} .= "$start";
    $prot_hash{$prot}{stop} .= "$stop";
    $prot_hash{$prot}{size} .= $stop - $start;
  }
  $tmp_prot = $prot;
}
print "";
close IN or die "Can't close $gff: $!\n";
