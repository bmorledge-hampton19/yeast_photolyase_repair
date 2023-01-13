#!/usr/bin/perl

use strict;
use warnings;

my $header = <STDIN>;
chomp $header;

my @col = split /\t/, $header;

print "GID\t$col[0]\tNAME\tGWEIGHT";

for ( my $i = 1; $i < scalar @col; $i++ )
{
	print "\t$col[$i]";
}
print "\n";

print "EWEIGHT\t\t\t";
for ( my $i = 1; $i < scalar @col; $i++ )
{
        print "\t0";
}
print "\n";

my $genenum = 0;

while ( my $line = <STDIN> )
{
	chomp $line;
	my @fields = split /\t/, $line, -1;

	print "GENE" . $genenum . "X\t$fields[0]\t$fields[0]\t0";
	for ( my $i = 1; $i < scalar @fields; $i++ )
	{
        	print "\t$fields[$i]";
	}
	print "\n";
	$genenum++;
}
