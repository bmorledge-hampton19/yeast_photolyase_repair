#!/usr/bin/perl

use strict;
use warnings;

print STDERR "Please enter filename for data matrix\n";
my $filename = <STDIN>;
chomp $filename;

open( FILE, $filename ) || die "Couldn't open file: $filename\n";

my $upstream_offset = -500;
my $downstream_offset = 640;
my $totalbins = $downstream_offset - $upstream_offset + 1;
my @cpd_sum;
my @dipy_sum;
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
				$cpd_sum[$i - 2] += $fields[$i];
			}
		}
		elsif ( $fields[1] eq "Dipyrimidines" )
		{
                       for( my $i = 2; $i < scalar @fields; $i++ )
                       {
                                $dipy_sum[$i - 2] += $fields[$i];
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
print "TS";
for (my $i = $upstream_offset; $i <= $downstream_offset; $i++ )
{
	print "\t$i";
}
print "\nCPDs";
my $midway = (scalar @cpd_sum) / 2;
for (my $i = 0; $i < $midway; $i++)
{
	print "\t$cpd_sum[$i]";
}
print "\nDipyrimidines";
if ( $midway != scalar (@dipy_sum) / 2 )
{	
	die "Error NMP and Purine arrays are of different sizes!\n";
}
elsif ( $midway != $totalbins )
{
	die "Mismatch in number of bins in input matrix and sum_all file\n";
}
for (my $i = 0; $i < $midway; $i++)
{
        print "\t$dipy_sum[$i]";
}
print "\nNormalized CPDs";
for (my $i = 0; $i < $midway; $i++)
{
	my $avg = "";
	if ( $dipy_sum[$i] > 0 )
	{
		$avg = 1.0 * $cpd_sum[$i]/$dipy_sum[$i];
	}
        print "\t$avg";
}
print "\n\n";
print "NTS";
for (my $i = $upstream_offset; $i <= $downstream_offset; $i++ )
{
        print "\t$i";
}
print "\nCPDs";
for (my $i = $midway; $i < scalar @cpd_sum; $i++)
{
        print "\t$cpd_sum[$i]";
}
print "\nDipyrimidines";
for (my $i = $midway; $i < scalar @dipy_sum; $i++)
{
        print "\t$dipy_sum[$i]";
}
print "\nNormalized CPDs";
for (my $i = $midway; $i < scalar @cpd_sum; $i++)
{
	my $avg = "";
	if ( $dipy_sum[$i] > 0 )
	{
        	$avg = 1.0 * $cpd_sum[$i]/$dipy_sum[$i];
	}
        print "\t$avg";
}
print "\n\n";
print "Strand Avg";
for (my $i = $upstream_offset; $i <= $downstream_offset; $i++ )
{
        print "\t$i";
}
print "\nCPDs";
for (my $i = 0; $i < $midway; $i++)
{
	my $bothstrand_cpd = $cpd_sum[$i] + $cpd_sum[($midway + $i)];
        print "\t$bothstrand_cpd";
}
print "\nDipyrimidines";
for (my $i = 0; $i < $midway; $i++)
{
	my $bothstrand_dipys = $dipy_sum[$i] + $dipy_sum[($midway + $i)];
        print "\t$bothstrand_dipys";
}
print "\nNormalized CPDs";
for (my $i = 0; $i < $midway; $i++)
{
	my $bothcpd = $cpd_sum[$i] + $cpd_sum[($midway + $i)];
	my $bothdipy = $dipy_sum[$i] + $dipy_sum[($midway + $i)];
	my $avg = "";
	if ( $bothdipy > 0 )
	{
        	$avg = 1.0 * $bothcpd/$bothdipy;
	}
        print "\t$avg";
}

print "\n";
