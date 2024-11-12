#!/usr/bin/perl

use strict;
use warnings;

print "Enter name of cluster file to center\n";
my $clusterfile = <STDIN>;
chomp $clusterfile;

my $clusterout = $clusterfile;

print "What is the value for centering (fraction remaining value)?\n";
my $center = <STDIN>;
chomp $center;

$clusterout =~ s/\.txt/_center${center}\.txt/ || die "file not .txt file\n";
open ( CLUSTER, "$clusterfile" ) || die "Couldn't open file $clusterfile\n";

open ( OUT, ">$clusterout" ) || die "Couldn't open file $clusterout\n";

my $extraheader = <CLUSTER>;
print "1st line: $extraheader\n";
my $header = <CLUSTER>;
print OUT $header;

my @data = ();
my @centeredata = ();
while ( my $line = <CLUSTER> )
{
	chomp $line;
	
	# to keep empty last tab, see https://stackoverflow.com/questions/3711649/perl-split-with-empty-text-before-after-delimiters
	my @fields = split /\t/, $line, -1;

	print OUT $fields[0];
	for ( my $i = 1; $i < scalar @fields; $i++ )
	{
		my $val = "";
		if ( $fields[$i] ne "" && $fields[$i] != 0 )
		{
			$val = $fields[$i];
			push @data, $val;

			# center
			$val -= $center;
			push @centeredata, $val;
		}
		print OUT "\t$val";
	}
	print OUT "\n";
}

# calculate median, based off of code from: https://stackoverflow.com/questions/5119034/using-perl-to-find-median-mode-standard-deviation

my $median;
my $midindex = int ( (scalar @data) / 2 );
my @sorted = sort {$a <=> $b} @data;
if ( (scalar @data) % 2 == 1 )
{
	$median = 1.0 * $sorted[$midindex];
}
else
{
	$median = ( $sorted[$midindex] + $sorted[$midindex - 1])/ 2.0;
}


my $cenmedian;
$midindex = int ( (scalar @centeredata) / 2 );
@sorted = sort {$a <=> $b} @centeredata;
if ( (scalar @centeredata) % 2 == 1 )
{
        $cenmedian = 1.0 * $sorted[$midindex];
}
else
{
        $cenmedian = ( $sorted[$midindex] + $sorted[$midindex - 1])/ 2.0;
}

print STDERR "Original median is $median\nCentered median is $cenmedian\n";

