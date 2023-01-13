#!/usr/bin/perl

use strict;
use warnings;

use lib '../';

use NucCoord;
use CPDNucReads;

print STDERR "Enter filename for output .txt file\n";
my $outfile = <STDIN>;

chomp $outfile;

print STDERR "Enter window size in bp\n";
my $nucwindow = <STDIN>;
chomp $nucwindow;
my $bedwindow = $nucwindow;

$outfile =~ s/\.txt/_${nucwindow}bp\.txt/ || die "output file must be a .txt file\n";

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

print STDERR "Loading nucleosome coordinates\n";
my $nucs = NucCoord->new();

my %chrdyad = $nucs->get_chr_dyad();
my $nuc_filename = $nucs->get_nuc_filename();
#print OUT header
print OUT "Nucleosome positioning data from: $nuc_filename\nSequencing data from file: @plusfiles\t@minusfiles";

my @pluscpdval;
my @pluscpdcount;
my @minuscpdval;
my @minuscpdcount;
my $nucount = 0;

for ( my $a = 0; $a < scalar @plusfiles; $a++ )
{
        print STDERR "Loading Probe Values\n";
        my $reads = CPDNucReads->new($plusfiles[$a], $minusfiles[$a]);

	foreach my $chr (sort keys %chrdyad)
	{
		print STDERR "Starting $chr\n";
		my %plusreads = $reads->get_plus_reads_for_chromosome($chr);
		my $num_plusreads = scalar keys %plusreads;
		my %minusreads = $reads->get_minus_reads_for_chromosome($chr);
		my $num_minusreads = scalar keys %minusreads;
		print STDERR "$chr reads: $num_plusreads plus reads and $num_minusreads minus reads\n";
		foreach my $dyad ( @{$chrdyad{$chr}} )
		{
			my $start = $dyad - $nucwindow;
			my $end = $dyad + $nucwindow;
			for ( my $i = (-1 * $nucwindow); $i <= $nucwindow; $i++ )
			{	
				my $pos = $dyad + $i;
				my $j = $i + $nucwindow;
				if ( exists $plusreads{$pos} )
				{
					$pluscpdval[$j] += $plusreads{$pos};
					if ( $a == 0 )
					{
						$pluscpdcount[$j]++;
					}
				}
				
	                        if ( exists $minusreads{$pos} )
	                        {
	                                $minuscpdval[$j] += $minusreads{$pos};
					if ( $a == 0 )
					{
	                                	$minuscpdcount[$j]++;
					}
	                        }
			}
			if ( $a == 0 )
			{
                                my $bedstart = $dyad - $bedwindow;
                                my $bedend = $dyad + $bedwindow;

                                # make start coord 0-based
                                $bedstart--;
                        	print BED "$chr\t$bedstart\t$bedend\t${chr}:$dyad\t.\t+\n";	
				$nucount++;
			}
		}
	
	}
	undef $reads; 
}

if ( scalar @pluscpdval != scalar @minuscpdval )
{
        die "mismatched array counts\n";
}
print OUT "\tCount of nucleosomes: $nucount\n";
print OUT "Position relative to dyad:";

for ( my $i = (-1 * $nucwindow); $i <= $nucwindow; $i++ )
{
        print OUT "\t$i";
}
print OUT "\n";

print OUT "Plus strand";
for (my $i = (-1 * $nucwindow); $i <= $nucwindow; $i++ )
{
	my $j = $i + $nucwindow;
	my $mean = 1.0 * $pluscpdval[$j] / $pluscpdcount[$j];
	print OUT "\t$mean";
}
print OUT "\nMinus strand";
for (my $i = (-1 * $nucwindow); $i <= $nucwindow; $i++ )
{
	my $j = $i + $nucwindow;
        my $mean = 1.0 * $minuscpdval[$j] / $minuscpdcount[$j];
        print OUT "\t$mean";
}
print OUT "\n";
print OUT "\n5\'->3\'Minus strand";
for (my $i = $nucwindow; $i >= (-1 * $nucwindow); $i-- )
{
        my $j = $i + $nucwindow;
        my $mean = 1.0 * $minuscpdval[$j] / $minuscpdcount[$j];
        print OUT "\t$mean";
}

print OUT "\n\n";
print OUT "Plus strand\n";
print OUT "UVPPs:";
for (my $i = (-1 * $nucwindow); $i <= $nucwindow; $i++ )
{
        my $j = $i + $nucwindow;
        print OUT "\t$pluscpdval[$j]";
}
print OUT "\nPotential lesion sites:";
for (my $i = (-1 * $nucwindow); $i <= $nucwindow; $i++ )
{
        my $j = $i + $nucwindow;
        print OUT "\t$pluscpdcount[$j]";
}
print OUT "\nMinus strand\n";
print OUT "UVPPs:";
for (my $i = (-1 * $nucwindow); $i <= $nucwindow; $i++ )
{
        my $j = $i + $nucwindow;
        print OUT "\t$minuscpdval[$j]";
}
print OUT "\nPotential lesion sites:";
for (my $i = (-1 * $nucwindow); $i <= $nucwindow; $i++ )
{
        my $j = $i + $nucwindow;
        print OUT "\t$minuscpdcount[$j]";
} 
print OUT "\n";
print OUT "\n5\'->3\'Minus strand\n";
print OUT "UVPPs:";
for (my $i = $nucwindow; $i >= (-1 * $nucwindow); $i-- )
{       
        my $j = $i + $nucwindow;
        print OUT "\t$minuscpdval[$j]";
}
print OUT "\nPotential lesion sites:";
for (my $i = $nucwindow; $i >= (-1 * $nucwindow); $i-- )
{       
        my $j = $i + $nucwindow;
        print OUT "\t$minuscpdcount[$j]";
}

print OUT "\n\n";
print OUT "Average of both strands:";
print OUT "\nAverage";
for (my $i = (-1 * $nucwindow); $i <= $nucwindow; $i++ )
{
        my $j = $i + $nucwindow;
	my $total_val = $pluscpdval[$j] + $minuscpdval[$j];
	my $total_count = $pluscpdcount[$j] + $minuscpdcount[$j];
        my $mean = 1.0 * $total_val / $total_count;
        print OUT "\t$mean";
}
