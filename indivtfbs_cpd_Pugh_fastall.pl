#!/usr/bin/perl

use strict;
use warnings;

use lib '../';
use TfbsYeastCoord;
use CPDNucReads;

print STDERR "Enter filename for output .txt file\n";
my $outfile = <STDIN>;
chomp $outfile;

my $bed_file = $outfile;
$bed_file =~ s/\.txt/\.bed/ || die "output file must be a .txt file\n";
open (BED, ">$bed_file") || die "couldn't open bed file\n";

open ( OUT, ">$outfile" ) || die "couldn't open output file\n";

my $window = 100;

# adsk for multiple wig filenames to analyze
my @plusfiles;
my @minusfiles;

my $counter = 1;
while ( 1 )
{
        print STDERR "Enter filename of plus strand reads for replicate $counter\n";
        my $plusfile = <STDIN>;
        chomp($plusfile);

        if ( $plusfile =~ /wig/ )
        {
                push @plusfiles, $plusfile;
        }
        else
        {
                last;
        }

        print STDERR "Enter filename of minus strand reads for replicate $counter\n";
        my $minusfile = <STDIN>;
        chomp($minusfile);

        if ( $minusfile =~ /wig/ )
        {
                push @minusfiles, $minusfile;
        }
        else
        {
                last;
        }

        $counter++;
}

if ( scalar @plusfiles != scalar @minusfiles )
{
        die "Number of plus and minus files given don't match!\n";
}

# get offset file
my $offsetfile = "../yeast_tfmotif_offsets.txt";
open (OFFSET, "$offsetfile" ) || die "Couldn't open $offsetfile\n";

my %offset = ();
my $offheader = <OFFSET>;
while ( my $off = <OFFSET> )
{
        chomp $off;
        my @field = split /\t/, $off;
        $field[0] =~ s/\s//g;
        $offset{$field[0]} = $field[1];
}
print STDERR "Loading TF coordinates\n";
my $tfsites = TfbsYeastCoord->new();

my %tf = $tfsites->get_tf_boundaries();
my $tf_filename = $tfsites->get_tf_filename();

#print header
print OUT "TF binding site data from: $tf_filename\nSequencing data from file: @plusfiles\t@minusfiles\n";

my %peakcpds;
my %restcpds;

my $peakstart = -4;
my $peakend = 4;

my %numsites;

for ( my $a = 0; $a < scalar @plusfiles; $a++ )
{
        print STDERR "Loading Probe Values\n";
        my $reads = CPDNucReads->new($plusfiles[$a], $minusfiles[$a]);
        foreach my $chr (sort keys %tf)
        {
                print STDERR "Starting $chr\n";
                my %plusreads = $reads->get_plus_reads_for_chromosome($chr);
                my $num_plusreads = scalar keys %plusreads;
                my %minusreads = $reads->get_minus_reads_for_chromosome($chr);
                my $num_minusreads = scalar keys %minusreads;
                print STDERR "$chr reads: $num_plusreads plus reads and $num_minusreads minus reads\n";
                my @tfpos = @{$tf{$chr}};
                for ( my $i = 0; $i < scalar @tfpos; $i++ )
                {
	                my $tfstart = $tfpos[$i]->[0];
	                my $tfend = $tfpos[$i]->[1];
                        my $tfstrand = $tfpos[$i]->[2];
			my $tfname = $tfpos[$i]->[3];
			if ( $a == 0 )
			{	
				$numsites{$tfname}++;
			}
	
			# increase start position by 1, since wig files (data files) are one-based, and bed files (position files) are zero-based (only start, not end coordinate)
			$tfstart++;

			my $tfmidpoint = ( $tfstart + $tfend ) / 2.0;
	
			if ( $tfstrand eq "+" )
			{		
				if ( $tfmidpoint - int($tfmidpoint) != 0.0 )
				{
					# if fractional value, round up for both strands
					$tfmidpoint = int ( $tfmidpoint + 1);
				}
				
                                # correct for offset
				if ( exists $offset{$tfname} )
				{
                                	$tfmidpoint += $offset{$tfname};
				}
			}	
			elsif ( $tfstrand eq "-" )
			{
                                if ( $tfmidpoint - int($tfmidpoint) != 0.0 )
                                {
                                        # if fractional value, round up for both strands
                                        $tfmidpoint = int ( $tfmidpoint );
                                }

                                # correct for offset
                                if ( exists $offset{$tfname} )
                                {
                                        $tfmidpoint -= $offset{$tfname};
                                }
			}	
			else
			{
				die "Weird strand: $tfstrand\n";
			}			
			# expand size of window
			my $winstart = $tfmidpoint - $window;
			my $winend = $tfmidpoint + $window;

			if ( $a == 0 )
			{
				my $bedstart = $winstart - 1;
				# set all strands as + /default
		                print BED "$chr\t$bedstart\t$winend\t$tfname\t.\t$tfstrand\n";
			}
	
			for ( my $j = $winstart; $j <= $winend; $j++ )
			{	
				my $relpos = $j - $tfmidpoint;
	
				if ( exists $plusreads{$j} )
				{
					if ( $relpos >= $peakstart && $relpos <= $peakend )
					{
						$peakcpds{$tfname} += $plusreads{$j};
					}
					else
					{
						$restcpds{$tfname} += $plusreads{$j};
					}
				}
                                if ( exists $minusreads{$j} )
                                {
                                        if ( $relpos >= $peakstart && $relpos <= $peakend )
                                        {
                                                $peakcpds{$tfname} += $minusreads{$j};
                                        }
                                        else
                                        {
                                                $restcpds{$tfname} += $minusreads{$j};
                                        }
                                }
			}
		}
	}
}


# print header
print OUT "TF name\tNum. 6-4PPs at TF midpoint [$peakstart to $peakend]\tNum. 6-4PPs elsewhere\tNum. of TFBS\tPeak 6-4PP Density\tFlank 6-4PP Density\tPeak-to-Flank Ratio\tNorm. Ratio\n";
my $peakwidth = $peakend - $peakstart + 1;
my $windowidth = ( 2 * $window ) + 1;
my $restwidth = $windowidth - $peakwidth;
foreach my $tf (sort keys %numsites )
{
	my $sites = $numsites{$tf};
	my $peak = 0;
	my $rest = 0;
	if ( exists $peakcpds{$tf} )
	{
		$peak = $peakcpds{$tf};
	}
	if ( exists $restcpds{$tf} )
	{
		$rest = $restcpds{$tf};
	}
	my $ratio = 0;
	my $normratio = 0;
	if ( $rest != 0 )
	{
		$ratio = 1.0 * $peak / $rest;
		$normratio = $ratio * $restwidth / $peakwidth;
	}

	my $peakdensity = 1.0 * $peak / ( $sites * $peakwidth );
	my $restdensity = 1.0 * $rest / ( $sites * $restwidth );

	print OUT "$tf\t$peak\t$rest\t$sites\t$peakdensity\t$restdensity\t$ratio\t$normratio\n";
} 

=pod
my $fa_file = $bed_file;
$fa_file =~ s/\.bed/\.fa/;
system ("fastaFromBed -s -name -fi ../saccer3_genome.fa -bed $bed_file -fo $fa_file" );
=cut
