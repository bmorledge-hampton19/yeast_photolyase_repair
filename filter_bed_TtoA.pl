use strict;
use warnings;

while ( my $line = <STDIN> )
{
	chomp $line;
	my @fields = split /\t/, $line;
	my $orig = substr $fields[3], 1, 1;
	my $sub = "${orig}>$fields[4]";
	if ( $sub eq "T>A" )
	{
		print "$line\n";
	}
}
