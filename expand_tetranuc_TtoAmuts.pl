#!/usr/bin/perl

use strict;
use warnings;

while ( my $line = <STDIN> )
{
	chomp $line;
	my @fields = split /\t/, $line;
	if ( $fields[3] =~ /[ATGC][ATGC][ATGC]/ )
	{
		$fields[1] -= 1;
		$fields[2] += 1;
		my $st = $fields[5];
		# get tetranuc, 3' end of T>A mutatio
		if ( $st eq "+" )
		{
			$fields[2]++;
		}
		elsif ( $st eq "-" )
		{
			$fields[1]--;
		}
		else
		{
			die "Weird strand!\n";
		}
		
		print "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]\n";
	}
	else
	{
		die "error with line: $line\n";
	}

}
