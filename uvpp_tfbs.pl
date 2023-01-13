#!/usr/bin/perl

use strict;
use warnings;

use lib '../';
use TFCoord;
use CPDNucReads;

print STDERR "Enter filename for output .txt file\n";
my $outfile = <STDIN>;

chomp $outfile;

print STDERR "Enter window size in bp\n";
my $window = <STDIN>;
chomp $window;
my $bedwindow = $window;

$outfile =~ s/\.txt/_${window}bp\.txt/ || die "output file must be a .txt file\n";

#$window -= 0.5; # set to half integers for inbetween cpd-seq data

my $bed_file = $outfile;
$bed_file =~ s/\.txt/\.bed/ || die "output file must be a .txt file\n";
open (BED, ">$bed_file") || die "couldn't open bed file\n";
open ( OUT, ">$outfile" ) || die "couldn't open output file\n";

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

print STDERR "Loading TF coordinates\n";
my $tfsites = TFCoord->new();

my %tf = $tfsites->get_tf_boundaries();
my $tf_filename = $tfsites->get_tf_filename();
#print OUT header
print OUT "TF binding site data from: $tf_filename\nSequencing data from file: @plusfiles\t@minusfiles\t";

my @pluscpdval;
my @pluscpdcount;
my @minuscpdval;
my @minuscpdcount;

my $tfbscount = 0;

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
			# increase start positions by 1, since wig files (data files) are one-based, and bed files (position files) are zero-based (only start, not end coordinate)
			$tfstart++;
	
			if ( $tfstrand eq "+" )
			{
				my $tfmidpoint = ( $tfstart + $tfend ) / 2.0;

                                if ( $tfmidpoint - int($tfmidpoint) != 0.0 )
                                {
                                        $tfmidpoint = int ( $tfmidpoint + 1);
                                }

                        	# expand size of window
	                        my $windowstart = $tfmidpoint - $window;
	                        my $windowend = $tfmidpoint + $window;
	                        if ( $a == 0 )
	                        {
	                                my $bedstart = $tfmidpoint - $bedwindow;
	                                my $bedend = $tfmidpoint + $bedwindow;
	
	                                # make start coord 0-based
	                                $bedstart--;
	                                print BED "$chr\t$bedstart\t$bedend\t.\t.\t$tfstrand\n";
					$tfbscount++;	
	                        }
				for ( my $j = $windowstart; $j <= $windowend; $j++ )
				{	
					my $relpos = $j - $windowstart;
					if ( exists $plusreads{$j} )
					{
						$pluscpdval[$relpos] += $plusreads{$j};
						if ( $a == 0 )
						{
							$pluscpdcount[$relpos]++;
						}
					}
				
	        	                if ( exists $minusreads{$j} )
	                	        {
	                        	        $minuscpdval[$relpos] += $minusreads{$j};
                                                if ( $a == 0 )
                                                {
	                                		$minuscpdcount[$relpos]++;
						}
		                        }
				}
			}
			elsif ( $tfstrand eq "-" )
			{
                                my $tfmidpoint = ( $tfstart + $tfend ) / 2.0;

                                if ( $tfmidpoint - int($tfmidpoint) != 0.0 )
                                {
                                        $tfmidpoint = int ( $tfmidpoint );
                                }

                                # expand size of window
                                my $windowstart = $tfmidpoint - $window;
                                my $windowend = $tfmidpoint + $window;
                                if ( $a == 0 )
                                {       
                                        my $bedstart = $tfmidpoint - $bedwindow;
                                        my $bedend = $tfmidpoint + $bedwindow;

                                        # make start coord 0-based
                                        $bedstart--;

                                        print BED "$chr\t$bedstart\t$bedend\t.\t.\t$tfstrand\n";
					$tfbscount++;
                                }
	                        for ( my $j = $windowend; $j >= $windowstart; $j-- )
	                        {
	                                my $relpos = $windowend - $j;
	                                if ( exists $minusreads{$j} )
	                                {
	                                        $pluscpdval[$relpos] += $minusreads{$j};
                                                if ( $a == 0 )
                                                {
	                                        	$pluscpdcount[$relpos]++;
						}
	                                }
	
	                                if ( exists $plusreads{$j} )
	                                {
	                                        $minuscpdval[$relpos] += $plusreads{$j};
                                                if ( $a == 0 )
                                                {
	                                        	$minuscpdcount[$relpos]++;
						}
	                                }
	                        }
	
	
	
	
			}
		}
	
	}
	undef $reads;
}

if ( scalar @pluscpdval != scalar @minuscpdval )
{
	die "mismatched array counts\n";
}
print OUT "Total count of tfbs: $tfbscount\n";
print OUT "CPD read counts";
print OUT "\nPlus strand:";
for (my $i = 0; $i < scalar @pluscpdval; $i++ )
{
        print OUT "\t$pluscpdval[$i]"; 
}
print OUT "\nMinus strand:";
for (my $i = 0; $i < scalar @minuscpdval; $i++ )
{
        print OUT "\t$minuscpdval[$i]"; 
}       
print OUT "\n\nTotal dipyrimidine counts";
print OUT "\nPlus strand:";
for (my $i = 0; $i < scalar @pluscpdcount; $i++ )
{
        print OUT "\t$pluscpdcount[$i]"; 
}       
print OUT "\nMinus strand:";
for (my $i = 0; $i < scalar @minuscpdcount; $i++ )
{
        print OUT "\t$minuscpdcount[$i]";
}
print OUT "\n\nAverage normalized CPDs:";
print OUT "\nPlus strand:";
for (my $i = 0; $i < scalar @pluscpdval; $i++ )
{
	if ( $pluscpdcount[$i] != 0 )
	{
        	my $mean = 1.0 * $pluscpdval[$i] / $pluscpdcount[$i];
        	print OUT "\t$mean";
	}
	else
	{
		print OUT "\t ";
	}
}
print OUT "\nMinus strand:";
for (my $i = 0; $i < scalar @minuscpdval; $i++ )
{
	if ( $minuscpdcount[$i] != 0 )
	{
        	my $mean = 1.0 * $minuscpdval[$i] / $minuscpdcount[$i];
        	print OUT "\t$mean";
	}
	else
	{
		print OUT "\t ";
	}
}
print OUT "\n\nPosition Relative to Motif Midpoint";
for (my $i = 0; $i < scalar @pluscpdcount; $i++ )
{
        my $pos = $i - $window;
        print OUT "\t$pos";
}
#both strands
print OUT "\n\nBoth strands\n";
print OUT "Total CPD count";
for (my $i = 0; $i < scalar @pluscpdval; $i++ )
{
	my $total = 0;
	$total = $pluscpdval[$i] + $minuscpdval[$i];
        print OUT "\t$total";
}
print OUT "\nTotal Dipyrimidine count";
for (my $i = 0; $i < scalar @pluscpdcount; $i++ )
{
        my $total = 0;
        $total = $pluscpdcount[$i] + $minuscpdcount[$i];
        print OUT "\t$total";
}
print OUT "\n\nNormalized CPDs";
for (my $i = 0; $i < scalar @pluscpdval; $i++ )
{
        my $totalcpds = 0;
        $totalcpds = $pluscpdval[$i] + $minuscpdval[$i];
	my $totaldipy = 0;
        $totaldipy = $pluscpdcount[$i] + $minuscpdcount[$i];
	if ( $totaldipy != 0 )
	{
		my $mean = 1.0 * $totalcpds / $totaldipy;
		print OUT "\t$mean";
	}
	else
	{
		print OUT "\t ";
	}
}
print OUT "\n\nPosition Relative to Motif Midpoint";
for (my $i = 0; $i < scalar @pluscpdcount; $i++ )
{
        my $pos = $i - $window;
        print OUT "\t$pos";
}

# create fasta file from bed
my $fa_file = $bed_file;
$fa_file =~ s/\.bed/\.fa/;
#system ("fastaFromBed -s -name -fi ../saccer3_genome.fa -bed $bed_file -fo $fa_file" );
