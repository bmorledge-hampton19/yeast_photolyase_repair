#!/usr/bin/perl

use strict;
use warnings;

package CPDNucReads;

my %chrlookup = ( "chr1" => "chrI", "chr2" => "chrII", "chr3" => "chrIII", "chr4" => "chrIV", "chr5" => "chrV", "chr6" => "chrVI", "chr7" => "chrVII", "chr8" => "chrVIII", "chr9" => "chrIX", "chr10" => "chrX", "chr11" => "chrXI", "chr12" => "chrXII", "chr13" => "chrXIII", "chr14" => "chrXIV", "chr15" => "chrXV", "chr16" => "chrXVI", "chrM" => "chrM");

sub new
{
	my ($class, $plusfile, $minusfile) = @_;

	my $self = bless {}, $class;

	my $inChr_flag = 0;
	my $sequence_num = 0;
	my $chr_num;
	my %plusreads;
	my %minusreads;

	open( PLUS, $plusfile) || die "couldn't open file: $plusfile\n";

	while( my $line = <PLUS> )
	{
	        chomp $line;
	        $line =~ s/^M//g;
		if ( $line eq "" )
		{
			next;
		}
	        elsif( $line =~ /chrom=(chr[IXVM]+)/ )
	        {
	                $chr_num = $1;
	        }
	        else
	        {
			my @fields = split /\t/, $line;

			# hash of hashes
	                $plusreads{$chr_num}{$fields[0]} = $fields[1];
		
			#print "field0: $fields[0]\t field1: $fields[1]\n";
	        }
	
	}
	close( PLUS );

	$self->{'plusreads'} = \%plusreads;

	open( MINUS, $minusfile) || die "couldn't open file: $minusfile\n";
	$chr_num = "";
        while( my $line = <MINUS> )
        {
                chomp $line;
                $line =~ s/^M//g;
                if ( $line eq "" )
                {
                        next;
                }
                elsif( $line =~ /chrom=(chr[IXVM]+)/ )
                {
                        $chr_num = $1;
                }
                else
                {
                        my @fields = split /\t/, $line;

                        $minusreads{$chr_num}{$fields[0]} = $fields[1];

                        #print "field0: $fields[0]\t field1: $fields[1]\n";
                }

        }
        close( MINUS );

        $self->{'minusreads'} = \%minusreads;

	return $self;
}

sub get_minus_reads_for_chromosome
{
	my( $self, $chr_num ) = @_;

 		
	return %{$self->{'minusreads'}{$chr_num}};
}

sub get_plus_reads_for_chromosome
{
        my( $self, $chr_num ) = @_;

        return %{$self->{'plusreads'}{$chr_num}};
}
1;
