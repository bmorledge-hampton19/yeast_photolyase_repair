#!/usr/bin/perl

use strict;
use warnings;

print STDERR "What is the name of the .bed file containing tfbs window coordinates?\n";
my $origfile = <STDIN>;

my $nucwindow = 100;
my $totalwindow = ( 2 * $nucwindow ) + 1;
chomp $origfile;
open( ORIG, $origfile ) || die "Couldn't open bed file: $origfile\n";

my $bedfile = $origfile;
$bedfile =~ s/\.bed/_expand\.bed/ || die "File is not bed file!\n";
open ( BED, ">$bedfile" ) || die "Couldn't open file $bedfile\n";

# process, bed file, expanding it by 1 nucleotide on each side to get trinucs from -73 to +73
while ( my $line = <ORIG> )
{
	chomp $line;
	if ( $line eq "" )
	{
		next;
	}
	my @fields = split /\t/, $line;
	my $start = $fields[1] - 1;
	my $end = $fields[2] + 1;
	print BED "$fields[0]\t$start\t$end\t$fields[3]\t$fields[4]\t$fields[5]\n";
}

# create fasta file from bed
my $fa_file = $bedfile;
$fa_file =~ s/\.bed/\.fa/;
system ("fastaFromBed -s -name -fi ../saccer3_genome.fa -bed $bedfile -fo $fa_file" );

# open trinuc_frequency file
print STDERR "What is the name of the background file?\n";
#my $trifile = "rad16_UV_bothruns_Muts_allSNVs_sorted_notTtoA_trinucmutation_strand_background.txt";
my $trifile = <STDIN>;
chomp $trifile;
open (TRI, $trifile ) || die "Couldn't open file: $trifile\n";

my %trinucs;
# get rid of headers;
my $head1 = <TRI>;
my $head2 = <TRI>;
my $tricount = 0;
while ( my $line = <TRI> )
{
	chomp $line;
	my @fields = split /\t/, $line;
	if ( $fields[0] =~ /^([ATGC]{3})$/ )
	{
		my $context = $1;
		$trinucs{$context} = $fields[3];
		$tricount++;
	}
	else
	{
		die "Misformatted line: $line\n";
	}
}

if ( $tricount != 32 )
{
	die "Total number of trinucleotide contexts is $tricount\n";
}

# process sequence in fasta file
open(FASTA, $fa_file ) || die "Couldn't open fasta file\n";
my @mutfraction;
my $totalfraction = 0;
my @nucount;
my $oddtrinucount = 0;
while ( my $line = <FASTA> )
{
        chomp $line;
	if ( $line =~ /^>/ )
	{
		next;	# header of fasta
	}

	if ( length($line) != ( $totalwindow + 2 ) )
	{
		die "FASTA line of wrong length: $line\n";
	}
	
	for ( my $i = 0; $i < (length($line) - 2); $i++ )
	{
		my $triseq = substr( $line, $i, 3 );
		my $mid = substr( $triseq, 1, 1 );
		if ( $mid eq "G" || $mid eq "A" )
		{
			$triseq = reverse $triseq;
			$triseq =~ tr/ACGT/TGCA/;
		}

		if ( exists $trinucs{$triseq} )
		{
			$mutfraction[$i] += $trinucs{$triseq};	
			$totalfraction += $trinucs{$triseq};
			$nucount[$i]++;
		}
		else
		{
			die "No data for trinuc: $triseq\n";
			$oddtrinucount++;
		}
	}
}

my $outfile = $origfile;
$outfile =~ s/\.bed/_background\.txt/;
open (OUT, ">$outfile" ) || die "Couldn't open OUT file: $outfile\n";

print STDERR "Number of odd (i.e., N-containing) trinucleotides: $oddtrinucount\n";
print OUT "Sequence coordinates from $origfile\n";
if ( scalar @mutfraction == $totalwindow )
{
	print OUT "Position Relative to Motif Midpoint (bp)";
	for ( my $i = -1 * $nucwindow; $i <= $nucwindow; $i++ )
	{
		print OUT "\t$i";
	}
	print OUT "\n";
}
else
{
	die "Error with size of mutfraction array\n";
}

=pod
print OUT "Expected Fraction of Mutations in Window";
for ( my $i = 0; $i < scalar @mutfraction; $i++ )
{
	my $fraction = 1.0 * $mutfraction[$i] / $totalfraction;	
	print OUT "\t$fraction";
}
print OUT "\n";
=cut

print OUT "Expected Mutations";
for ( my $i = 0; $i < scalar @mutfraction; $i++ )
{
        my $mutation = $mutfraction[$i];
        print OUT "\t$mutation";
}
print OUT "\n";

print OUT "Total Nucleotides";
for ( my $i = 0; $i < scalar @nucount; $i++ )
{
        print OUT "\t$nucount[$i]";
}
print OUT "\n\n";

print OUT "Mutation Density";
for ( my $i = 0; $i < scalar @mutfraction; $i++ )
{
        my $mutdensity = 1.0 * $mutfraction[$i] / ( $nucount[$i] );
        print OUT "\t$mutdensity";
}
print OUT "\n";

