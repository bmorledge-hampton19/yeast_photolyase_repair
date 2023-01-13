#!/usr/bin/perl

use strict;
use warnings;

print STDERR "Opening FASTA file\n";
my $fa_file = "../saccer3_genome.fa";
# process sequence in fasta file
open(FASTA, $fa_file ) || die "Couldn't open fasta file\n";
print STDERR "Finished opening FASTA file\n";
my $totalnuc = 0;
my %tetranucs;
my $chr = "";
my $flag = 0;
my $seq = "";
while ( my $line = <FASTA> )
{
        chomp $line;
	if ( $line =~ /^>([a-zA-Z0-9_]+)/ )
	{
		# count tetranucs in sequence from previous chromosome
		for ( my $i = 0; $i < (length($seq) - 3); $i++ )
        	{
                	my $tetraseq = substr( $seq, $i, 4 );
                	$tetranucs{$tetraseq}++;
			$totalnuc++;

			# reverse comp
			$tetraseq =~ tr/ACGT/TGCA/;
			$tetraseq = reverse $tetraseq;
                        $tetranucs{$tetraseq}++;
                        $totalnuc++;

        	}
		$chr = $1;
		$seq = "";
		if ( $chr eq "pUC19" || $chr eq "chrM" )
		{
			print STDERR "$chr sequence was skipped\n";
			$flag = 1;
		}
		elsif ( $chr =~ /_/ )
		{
			print STDERR "Skipping chromosome $chr\n";
			$flag = 1;
		}
		else
		{
			print STDERR "Starting chromosome $chr\n";
			$flag = 0;
		}
		next;	# header of fasta
	}

	if ( $flag == 1 )
	{
		next;
	}
	else
	{
		$seq .= $line;
	}
}

if ( $flag == 0 )
{

        for ( my $i = 0; $i < (length($seq) - 3); $i++ )
        {
                my $tetraseq = substr( $seq, $i, 4 );
                $tetranucs{$tetraseq}++;
		$totalnuc++;

                        # reverse comp
                        $tetraseq =~ tr/ACGT/TGCA/;
                        $tetraseq = reverse $tetraseq;
                        $tetranucs{$tetraseq}++;
                        $totalnuc++;
        }
}
my $outfile = "saccer3_noChrM_tetranuc_background.txt";
open (OUT, ">$outfile" ) || die "Couldn't open OUT file: $outfile\n";

print OUT "Data from bed file: saccer3_genome.fa\n";
print OUT "Total nucleotide count: $totalnuc\n";

print OUT "Tetranucleotide context\tNucleotide count\tFraction of total nucleotides\n";
foreach my $tetra ( sort keys %tetranucs )
{
	my $fraction = 1.0 * $tetranucs{$tetra} / $totalnuc;
	print OUT "$tetra\t$tetranucs{$tetra}\t$fraction\n";
}
