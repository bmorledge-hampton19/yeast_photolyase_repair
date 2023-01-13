#!/usr/bin/perl

use strict;
use warnings;

print STDERR "Enter name of bed file:\n";
my $bedfile = <STDIN>;
chomp $bedfile;
open (BED, "$bedfile") || die "Could not open .bed file: $bedfile\n";

print STDERR "Enter name of matching dinucleotide sequence file:\n";
my $dinucfile = <STDIN>;
chomp $dinucfile;
open (DINUC, "$dinucfile") || die "Could not open .fa file: $dinucfile\n";

my $head = "";
my $dinucseq = "";
my $match_flag = 1;
while( my $line = <BED> )
{
	my @field = split /\t/, $line;
	if ( $match_flag )
	{	
		my $temp = <DINUC>;
		chomp $temp;
		$dinucseq = <DINUC>;
		chomp $dinucseq;
		if ( $temp =~ /\(/ )
		{
			my @fields = split /\(/, $temp;
			$head = $fields[0];
			$head =~ s/^>//;
		}
		else 
		{
			die "Misformatted line: $temp\n";
		}

	}
	if ( $head eq $field[3] )
	{
		$match_flag = 1;
		if ( $dinucseq =~ /[T][A]/ )
		{
			print $line;
		}
	}
	else
	{
		print STDERR "Error: sequence read names didn't match: $field[3] and $head\n";
		$match_flag = 0;
	}
}
