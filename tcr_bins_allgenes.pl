#!/usr/bin/perl

use strict;
use warnings;

print STDERR "Please enter filename for data matrix\n";
my $filename = <STDIN>;
chomp $filename;

open( FILE, $filename ) || die "Couldn't open file: $filename\n";
# number of bins to compute
my $numgenebins = 6;
my $numflankbins = 3;
my $interval = 167; # size of flanking interval bins
my $flank_offset = $numflankbins * $interval;
my $totalbins = $numgenebins + 2 * $numflankbins;

my @nmp_sum;
my @purine_sum;
while ( my $line = <FILE> )
{
	chomp $line;
	if ( $line =~ /^(Y[A-P][LR][0-9]{3}[CW]\-?[A-H]?)/ )
	{
		my $acc = $1;
		my @fields = split /\t/, $line;	

                if ( scalar @fields != ( 2 * $totalbins + 2 ) )
                {
                        die "Wrong number of bins for gene: $acc\n";
                }
		
		if ( $fields[1] eq "CPDs" )
		{
			for( my $i = 2; $i < scalar @fields; $i++ )
			{
				$nmp_sum[$i - 2] += $fields[$i];
			}
		}
		elsif ( $fields[1] eq "Dipyrimidines" )
		{
                       for( my $i = 2; $i < scalar @fields; $i++ )
                       {
                                $purine_sum[$i - 2] += $fields[$i];
                       }
                }
		else
		{
			die "Count type for acc $acc is $fields[1]\n";
		}
	}		
	else
	{
		print STDERR "No match for line: $line\n";
	}
}

print "From file: $filename\t for transcription frequency subsets\n";
# print results:

for (my $i = 0; $i < $numflankbins; $i++ )
{
        my $start = 0 - $flank_offset + ( $i * $interval);
        my $end = 0 - $flank_offset + (($i + 1) * $interval) - 1;
	print "\tTS Promoter ($start to $end)";
}
for ( my $i = 1; $i <= $numgenebins; $i++ )
{
	print "\tTS Coding bin $i";
}
for ( my $i = 0; $i < $numflankbins; $i++ )
{
	my $start = ( $i * $interval) + 1; 
        my $end = ($i + 1) * $interval;
	print "\tTS Terminator ($start to $end)";
}

print "\nCPDs";
my $midway = (scalar @nmp_sum) / 2;
for (my $i = 0; $i < $midway; $i++)
{
	print "\t$nmp_sum[$i]";
}
print "\nDipyrimidines";
if ( $midway != scalar (@purine_sum) / 2 )
{	
	die "Error NMP and Purine arrays are of different sizes!\n";
}
elsif ( $midway != ($numgenebins + 2 * $numflankbins ) )
{
	die "Mismatch in number of bins in input matrix and sum_all file\n";
}
for (my $i = 0; $i < $midway; $i++)
{
        print "\t$purine_sum[$i]";
}
print "\nNormalized CPDs";
for (my $i = 0; $i < $midway; $i++)
{
	my $avg = 1.0 * $nmp_sum[$i]/$purine_sum[$i];
        print "\t$avg";
}
print "\n\n";
for (my $i = 0; $i < $numflankbins; $i++ )
{
        my $start = 0 - $flank_offset + ( $i * $interval);
        my $end = 0 - $flank_offset + (($i + 1) * $interval) - 1;
        print "\tNTS Promoter ($start to $end)";
}
for ( my $i = 1; $i <= $numgenebins; $i++ )
{
        print "\tNTS Coding bin $i";
}
for ( my $i = 0; $i < $numflankbins; $i++ )
{
        my $start = ( $i * $interval) + 1; 
        my $end = ($i + 1) * $interval;
        print "\tNTS Terminator ($start to $end)";
}
print "\nCPDs";
for (my $i = $midway; $i < scalar @nmp_sum; $i++)
{
        print "\t$nmp_sum[$i]";
}
print "\nDipyrimidines";
for (my $i = $midway; $i < scalar @purine_sum; $i++)
{
        print "\t$purine_sum[$i]";
}
print "\nNormalized CPDs";
for (my $i = $midway; $i < scalar @nmp_sum; $i++)
{
        my $avg = 1.0 * $nmp_sum[$i]/$purine_sum[$i];
        print "\t$avg";
}
print "\n\n";
