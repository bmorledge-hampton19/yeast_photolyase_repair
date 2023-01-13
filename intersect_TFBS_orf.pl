#!/usr/bin/perl 
# intersect_TFBS.pl

use strict;
use warnings;

my $buffer = 100;
print STDERR "Enter filename of TFBS to intersect with DHS coordinates\n";
my $tfbs_file = <STDIN>;
chomp $tfbs_file;
open ( TFBS, $tfbs_file ) || die "Couldn't open file: $tfbs_file\n";

# open DHS coordinates
print STDERR "Enter filename of DHS coordinate file:\n";
my $dhs_file = <STDIN>;
chomp $dhs_file;
open (DHS, $dhs_file) || die "Couldn't open promoter file: $dhs_file\n";
print STDERR "processing DHS coordinates...\n";

# process promoters into lookup hash (1 => promoter, undef/not exist ==> not promoter)
# I am assuming these coordinates are 1 based
my %DHS_coord;
my $chr = "";
while ( my $line = <DHS> )
{
	chomp $line;
	my @fields = split /\t/, $line;
	if ( $fields[0] =~ /(chr[XIV]+)/ )
	{
		my $temp = $1;
		if ( $chr ne $temp )
		{
			print STDERR "Starting to process $temp\n";
			$chr = $temp;
		}

		push @{$DHS_coord{$chr}}, [$fields[1], $fields[2]];
	}
	else
	{
		#die "Misformatted line: $line\n";
	}
}
close (DHS);

# intersect TF coordinates with DHS --> any overlap means active DHS TF

$tfbs_file =~ s/\.bed//;

my $InDHS_file = $tfbs_file . "_ORF100bp.bed";

my $notInDHS_file = $tfbs_file . "_notORF100.bed";

open (DHS, ">$InDHS_file") || die "Couldn't open InDHS file: $InDHS_file\n";
open (NOTDHS, ">$notInDHS_file") || die "Couldn't open NOTInDHS file: $notInDHS_file\n";

print STDERR "Starting processing TFBS coordinates...\n";
$chr = "";
my %DHS_lookup;
while ( my $line = <TFBS>)
{
	chomp $line;
	my @fields = split /\t/, $line;
        if ( $fields[0] =~ /(chr[XIV]+)/ )
        {
		my $temp = $1;
		if ( $temp ne $chr )
		{
                	print STDERR "Starting to process $temp\n";
			$chr = $temp;
			# empty hash
			%DHS_lookup = ();
			foreach my $dhs ( @{$DHS_coord{$chr}} )
			{
				my $start = $dhs->[0];
				# start coordinate is 0-based, so convert to 1-based for consistency
				$start++;
				$start += $buffer;
				my $end = $dhs->[1];
				$end -= $buffer;
                		for ( my $i = $start; $i <= $end; $i++ )
                		{
                        		# Define this nucleotide as a promoter
                        		$DHS_lookup{$fields[0]}{$i} = 1;
                		}
			}
		}
                my $start = $fields[1] + 1; # to make it 1-based instead of 0-based for consistency
		my $end = $fields[2]; # already 1 based
		my $prox_flag = 0;
		for ( my $i = $start; $i <= $end; $i++ )
		{
			if ( exists $DHS_lookup{$chr}{$i} && $DHS_lookup{$chr}{$i} == 1 )
			{
				$prox_flag = 1;	
			}
		}

                if ( $prox_flag )
                {
                        print DHS "$line\n";
                }
                else
                {
                        print NOTDHS "$line\n";
                }
	}
	else 
	{
		#die "Misformatted line: $line\n";
	}
}
