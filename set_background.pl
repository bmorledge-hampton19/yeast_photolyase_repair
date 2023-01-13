#!/usr/bin/perl

use warnings;
use strict;

print STDERR "What is the name of the .wig file?\n";
my $wigfile = <STDIN>;
chomp $wigfile;

print STDERR "What is the name of the background file\n";
my $bkgdfile = <STDIN>;
chomp $bkgdfile;

open (BKGD, $bkgdfile ) || die "Couldn't open file $bkgdfile\n";
my $chr = "";
my %reads;
while ( my $line = <BKGD> )
{
	chomp $line;
	if ( $line =~ /(chr[IXVM]+)/ )
	{
		$chr = $1;
	}
	else
	{
		my @field = split /\t/, $line;
		$reads{$chr}{$field[0]} = $field[1];
	}

}
close ( BKGD );

open (WIG, $wigfile ) || die "Couldn't open file $wigfile\n";
$chr = "";
while ( my $line = <WIG> )
{
        chomp $line;
        if ( $line =~ /(chr[IXVM]+)/ )
        {
                $chr = $1;
        }
        else
        {
                my @field = split /\t/, $line;
                $reads{$chr}{$field[0]} = $field[1];
        }

}

foreach my $chrom (sort keys %reads )
{
        print "variableStep chrom=$chrom span=1\n";
	my %temp = %{$reads{$chrom}};
	foreach my $pos (sort { $a <=> $b } keys %temp )
	{
		print "$pos\t$temp{$pos}\n";
	}

} 
