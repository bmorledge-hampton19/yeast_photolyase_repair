#!/usr/bin/perl

use strict;
use warnings;

while ( my $line = <STDIN> )
{
	chomp $line;
	my @fields = split /\t/, $line;

	# expand to trinuc context
	my $start = $fields[1] + 96;
	my $end = $fields[2] - 96;

	if ($start<0)
	{ next;	}	
	print "$fields[0]\t$start\t$end\t$fields[3]\t$fields[4]\t$fields[5]\n";
}
