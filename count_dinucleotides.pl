#!/usr/bin/perl

use strict;
use warnings;

my %dinuc_count;

while ( my $line = <STDIN> )
{
	chomp $line;
	if ($line =~ /^>/ )
	{
		next;
	}
	
	$dinuc_count{$line}++;

}

foreach my $key (sort keys %dinuc_count)
{
	print "$key\t$dinuc_count{$key}\n";

}
