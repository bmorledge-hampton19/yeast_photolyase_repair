#!/usr/bin/perl

use strict;
use warnings;

use lib '../';
use GeneCoord;
use CPDReadValues;

# ask for probe filename to analyze
print STDERR "Enter filename of plus strand reads\n";
my $plusfile = <STDIN>;
chomp($plusfile);

print STDERR "Enter filename of minus strand reads\n";
my $minusfile = <STDIN>;
chomp($minusfile);

print STDERR "Loading Gene coordinates\n";
my $genes = GeneCoord->new();
print STDERR "Loading Probe Values\n";
my $reads = CPDReadValues->new($plusfile, $minusfile);

# location offsets
my $upstream_offset = -500;
my $downstream_offset = 640;

# START FROM HERE
my %chromosomes = $genes->get_chromosomes();
my %trxstart = $genes->get_tss();
my %trxend = $genes->get_tts();
my %strand = $genes->get_strand();

#print header
print "Data from file: $plusfile\t$minusfile\n";
print "Gene\tCOUNT_TYPE";

for (my $i = $upstream_offset; $i <= $downstream_offset; $i++ )
{
	print "\t$i (TS)";
}
for (my $i = $upstream_offset; $i <= $downstream_offset; $i++ )
{
        print "\t$i (NTS)";
}
print "\n";

foreach my $chr (sort keys %chromosomes)
{
	print STDERR "Starting $chr\n";
	my %plusreads = $reads->get_plus_reads_for_chromosome($chr);
	my $num_plusreads = scalar keys %plusreads;
	my %minusreads = $reads->get_minus_reads_for_chromosome($chr);
	my $num_minusreads = scalar keys %minusreads;
	print STDERR "$chr reads: $num_plusreads plus reads and $num_minusreads minus reads\n";
	foreach my $acc ( @{$chromosomes{$chr}} )
	{
		my $tss = $trxstart{$acc};
		my $tts = $trxend{$acc};
		
		my @cpd_ts = ();
		my @cpd_nts = ();
		my @dipy_ts = ();
		my @dipy_nts = ();
		# calculate read sums (CPD and DIPY bkgd) for gene
		for ( my $i = $upstream_offset; $i <= $downstream_offset; $i++)
		{
			my $pos;
			if ( $strand{$acc} eq "+" )
			{
				$pos = $tss + $i;
			}
			elsif ( $strand{$acc} eq "-" )
			{
				$pos = $tss - $i;
			}
			else
			{
				die "No strand information for gene: $acc\n";
			}
			my $plus_cpds = 0;
			my $plus_dipys = 0;
			my $minus_cpds = 0;
			my $minus_dipys = 0;
			if ( exists $plusreads{$pos} )
			{
				$plus_cpds += $plusreads{$pos};
				$plus_dipys++;
			}
			if ( exists $minusreads{$pos} )
			{
				$minus_cpds += $minusreads{$pos};
				$minus_dipys++;
			}
                        if ( $strand{$acc} eq "+" )
                        {
				push @cpd_nts, $plus_cpds;
				push @dipy_nts, $plus_dipys;

				push @cpd_ts, $minus_cpds;
				push @dipy_ts, $minus_dipys;
			}
			elsif ( $strand{$acc} eq "-" )
                        {
                                push @cpd_ts, $plus_cpds;
                                push @dipy_ts, $plus_dipys;

                                push @cpd_nts, $minus_cpds;
                                push @dipy_nts, $minus_dipys;
			}
		}

		if ( scalar @cpd_ts != scalar @cpd_nts || scalar @cpd_ts != scalar @dipy_ts || scalar @cpd_nts != scalar @dipy_nts )
		{
			die "Arrays are of different sizes!\n";
		}
  
		# print average probe values for acc
		print "$acc\tCPDs";

		foreach my $val (@cpd_ts)
		{
			print "\t$val";
		}
		foreach my $val (@cpd_nts)
		{
			print "\t$val";
		}
		print "\n";
		print "$acc\tDipyrimidines";
                foreach my $val (@dipy_ts)
                {
                        print "\t$val";
                }
                foreach my $val (@dipy_nts)
                {
                        print "\t$val";
                }
		print "\n";
	}

}
