#!/usr/bin/perl

use strict;
use warnings;

package NucCoord;

sub new
{
	my ($class) = @_;
	
	my $self = bless {}, $class;

	#my $nuc_file_name = "../H3K79me1_top10K.txt";
	#my $nuc_file_name = "../H3K79me1_bottom10K.txt";
	#my $nuc_file_name = "../H3K36me3_bottom10K.txt";
	#my $nuc_file_name = "../H3K36me3_top10K.txt";
	#my $nuc_file_name = "../H3K14Ac_top10K.txt";
	#my $nuc_file_name = "../H3K14Ac_bottom10K.txt";
	#my $nuc_file_name = "widomdyads_saccer3_1based_all.txt";
        my $nuc_file_name = "../widom_sac3_nucscore5_dyads.txt";
	#my $nuc_file_name = "../widom_sac3_nucscoreless1_dyads.txt";
	#my $nuc_file_name = "../Widom_dyads_saccer3_top2000_1based.txt";

	# open file with gene positions
	open( GENE, $nuc_file_name) || die "Couldn't open file\n";
	#my $header = <GENE>;
	
	my %chromdyad;
	
	while( my $line = <GENE> )
	{
		chomp($line);

		if ( $line =~ /(chr[XIVM]+):([0-9]+)-/ )
		{	
			push @{$chromdyad{$1}}, $2;
		}
	}

	close ( GENE );
	$self->{'chromdyad'} = \%chromdyad;
	$self->{'nuc_file_name'} = $nuc_file_name;
	return $self;	
}

sub get_chr_dyad
{
	my ($self) = @_;

	return %{$self->{'chromdyad'}};

}

sub get_nuc_filename
{
        my ($self) = @_;

        return $self->{'nuc_file_name'};
}

1;	
