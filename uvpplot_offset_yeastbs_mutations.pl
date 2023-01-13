#!/usr/bin/perl

use strict;
use warnings;

use lib '../';
use TfbsYeastCoord;
use CPDNucReads;

print STDERR "Enter filename for output .txt file\n";
my $outfile = <STDIN>;

chomp $outfile;

print STDERR "Enter window size in bp\n";
my $window = <STDIN>;
chomp $window;
my $bedwindow = $window;

print STDERR "Enter name of TFBS to analyze\n";
my $tfbsname = <STDIN>;
chomp $tfbsname;
$outfile =~ s/\.txt/_${window}bp_${tfbsname}\.txt/ || die "output file must be a .txt file\n";

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

# get offset file
my $offsetfile = "../yeast_tfmotif_offsets.txt";
open (OFFSET, "$offsetfile" ) || die "Couldn't open $offsetfile\n";

my %tfoffset = ();
my $offheader = <OFFSET>;
while ( my $off = <OFFSET> )
{
	chomp $off;
	my @field = split /\t/, $off;
	$field[0] =~ s/\s//g;
	$tfoffset{$field[0]} = $field[1];
}
print STDERR "Loading TF coordinates\n";
my $tfsites = TfbsYeastCoord->new();

my %tfname = $tfsites->get_tfname_boundaries();
my $tf_filename = $tfsites->get_tf_filename();
my %tf;
if ( exists $tfname{$tfbsname} )
{
	%tf = %{$tfname{$tfbsname}};
	
}
else 
{
	die "$tfbsname doesn't exist in file $tf_filename\n";
}
my $offset = 0;
if ( exists $tfoffset{$tfbsname} )
{
	$offset = $tfoffset{$tfbsname};
}
else
{	die "$tfbsname doesn't exist in file: $offsetfile\n";	}

#print OUT header
print OUT "TF binding site data for $tfbsname from: $tf_filename\nSequencing data from file: @plusfiles\t@minusfiles\t";

my %pluscpdval;
my %minuscpdval;

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

				# correct for offset
				$tfmidpoint += $offset;

                        	# expand size of window
	                        my $windowstart = $tfmidpoint - $window;
	                        my $windowend = $tfmidpoint + $window;
				my $tfbsid = "$chr:$tfmidpoint $tfstrand";
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
						$pluscpdval{$tfbsid}{$relpos} += $plusreads{$j};
					}
					else
					{
						$pluscpdval{$tfbsid}{$relpos} += 0;
					}
				
	        	                if ( exists $minusreads{$j} )
	                	        {
	                        	        $minuscpdval{$tfbsid}{$relpos} += $minusreads{$j};
		                        }
					else
                                        {
                                                $minuscpdval{$tfbsid}{$relpos} += 0;
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

                                # correct for offset
                                $tfmidpoint -= $offset;

                                # expand size of window
                                my $windowstart = $tfmidpoint - $window;
                                my $windowend = $tfmidpoint + $window;
                                my $tfbsid = "$chr:$tfmidpoint $tfstrand";

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
	                                        $pluscpdval{$tfbsid}{$relpos} += $minusreads{$j};
	                                }
                                        else
                                        {
                                                $pluscpdval{$tfbsid}{$relpos} += 0;
                                        }
	
	                                if ( exists $plusreads{$j} )
	                                {
	                                        $minuscpdval{$tfbsid}{$relpos} += $plusreads{$j};
	                                }
                                        else
                                        {
                                                $minuscpdval{$tfbsid}{$relpos} += 0;
                                        }
	                        }
	
	
	
	
			}
		}
	
	}
	undef $reads;
}

if ( scalar keys %pluscpdval != scalar keys %minuscpdval )
{
	die "mismatched array counts\n";
}
print OUT "\nTFBS position";
for (my $i = -1 * $window; $i <= $window; $i++ )
{
        print OUT "\t$i";
}
#both strands
my $corrsize = 2 * $window + 1;
print OUT "\n";
foreach my $id ( sort keys %pluscpdval )
{
	if ( scalar keys %{$pluscpdval{$id}} != $corrsize || scalar keys %{$minuscpdval{$id}} != $corrsize)
	{
		die "mismatched hashes!\n";
	}

	print OUT "$id";
	for (my $i = 0; $i < $corrsize; $i++ )
	{
		my $total = 0;
		$total = $pluscpdval{$id}{$i} + $minuscpdval{$id}{$i};
        	print OUT "\t$total";
	}

	print OUT "\n";
}
