#!/usr/bin/perl

use strict;
use warnings;

package TFCoord;

sub new
{
	my ($class) = @_;
	
	my $self = bless {}, $class;

	my $tf_file_name = "../Reb1_10MNase_80mM_len50_less10_expanded.bed";
	#my $tf_file_name = "../Reb1_10MNase_80mM_len50_greater10_expanded.bed";
	#my $tf_file_name = "gal4_chipexo_expanded.bed";
	#my $tf_file_name = "../Abf1_10MNase_80mM_len50_sites_expanded.bed";
	#my $tf_file_name = "Reb1_10MNase_80mM_len50_sites_expanded.bed";
	#my $tf_file_name = "../Rhee_SuppData1_TATA.bed";
	#my $tf_file_name = "../PughTFchipexo_motifs_mcm1.bed";
	#my $tf_file_name = "../Rhee_SuppData1_TATAformat_dipy.bed";
        #my $tf_file_name = "../Rhee_SuppData1_TATAformat_NOTdipy.bed";

	open( GENE, $tf_file_name) || die "Couldn't open file\n";
	
	my %tf;

	# get chromosome length
	my %chrlen = ();
	my $chrfile = "../chrom_len_saccer3.txt";
	open (CHROM, $chrfile ) || die "couldn't open file\n";
	while ( my $ch = <CHROM> )
	{
		chomp $ch;
		my @val = split /\t/, $ch;
		$chrlen{$val[0]} = $val[1];
	}
	
	my %excluded = ( "chrXII	458619	458662	-" => 1, "chrXII	467756	467799	-" => 1, "chrXII	490145	490188	-" => 1, "chrXII	490280	490327	-" => 1);
	while( my $line = <GENE> )
	{
		chomp($line);
		if ( $line eq "" )
		{	next;	}
		elsif ( exists $excluded{$line} && $excluded{$line} == 1 )
		{
			print STDERR "Excluded line: $line\n";
			next;
		}

		my @fields = split /\t/, $line;

		# remove tfbs near telomere
		if ( $fields[1] < 1000 || ( $chrlen{$fields[0]} - $fields[1] ) < 1000  )
		{
			print STDERR "Excluded tfbs: $line\n";
			next;
		}

		push @{$tf{$fields[0]}}, [$fields[1], $fields[2], $fields[3]];
		#print STDERR "$line\n";
	}

	close ( GENE );
	$self->{'tf'} = \%tf;
	$self->{'tf_file_name'} = $tf_file_name;
	return $self;	
}

sub get_tf_boundaries
{
	my ($self) = @_;

	return %{$self->{'tf'}};

}

sub get_tf_filename
{
        my ($self) = @_;

        return $self->{'tf_file_name'};
}

1;	
