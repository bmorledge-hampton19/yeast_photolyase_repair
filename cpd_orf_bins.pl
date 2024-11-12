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

# number of bins to compute
my $numgenebins = 6;
my $numflankbins = 3;

my %chromosomes = $genes->get_chromosomes();
my %trxstart = $genes->get_tss();
my %trxend = $genes->get_tts();
my %strand = $genes->get_strand();

my $interval = 167; # size of flanking interval bins
my $flank_offset = $numflankbins * $interval;

#print header
print "Data from file: $plusfile\t$minusfile\n";
print "Gene\tCOUNT_TYPE";

for (my $i = 0; $i < $numflankbins; $i++ )
{
        my $start = 0 - $flank_offset + ( $i * $interval);
        my $end = 0 - $flank_offset + (($i + 1) * $interval) - 1;
	print "\tTS Promoter ($start to $end)";
}
for ( my $i = 1; $i <= $numgenebins; $i++ )
{
	print "\tTS Coding bin $i";
}
for ( my $i = 0; $i < $numflankbins; $i++ )
{
	my $start = ( $i * $interval) + 1; 
        my $end = ($i + 1) * $interval;
	print "\tTS Terminator ($start to $end)";
}
for (my $i = 0; $i < $numflankbins; $i++ )
{
        my $start = 0 - $flank_offset + ( $i * $interval);
        my $end = 0 - $flank_offset + (($i + 1) * $interval) - 1;
        print "\tNTS Promoter ($start to $end)";
}
for ( my $i = 1; $i <= $numgenebins; $i++ )
{
        print "\tNTS Coding bin $i";
}
for ( my $i = 0; $i < $numflankbins; $i++ )
{
        my $start = ( $i * $interval) + 1; 
        my $end = ($i + 1) * $interval;
        print "\tNTS Terminator ($start to $end)";
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
		my $tts;
		if ( exists $trxend{$acc} )
		{
			$tts = $trxend{$acc};
		}
		else
		{
			next;
		}
		my $genelength = abs( $tss - $tts );
		my $genebinspan = 1.0 * $genelength / $numgenebins;

		#print "bin span: $genebinspan\n";
		
		my @cpd_ts = ();
		my @cpd_nts = ();
		my @dipy_ts = ();
		my @dipy_nts = ();
		if ( $strand{$acc} eq "+" )
		{
			# calculate read sums (NMP and Purine bkgd) for promoter
			for ( my $i = 0; $i < $numflankbins; $i++)
			{
				my $start = $tss - $flank_offset + ( $i * $interval);
				my $end = $tss - $flank_offset + (($i + 1) * $interval) - 1; # -1 because we don't include end nuc
				my $plus_cpds = 0;
				my $plus_dipys = 0;
				my $minus_cpds = 0;
				my $minus_dipys = 0;
				for ( my $j = $start; $j <= $end; $j++ )
				{
					if ( exists $plusreads{$j} )
					{
						$plus_cpds += $plusreads{$j};
						$plus_dipys++;
					}
					if ( exists $minusreads{$j} )
					{
						$minus_cpds += $minusreads{$j};
						$minus_dipys++;
					}
				}
				push @cpd_nts, $plus_cpds;
				push @dipy_nts, $plus_dipys;

				push @cpd_ts, $minus_cpds;
				push @dipy_ts, $minus_dipys;
			}
				
			# calculate read sums for coding region
			for (my $i = 0; $i < $numgenebins; $i++ )
			{
				my $start = $i * $genebinspan + $tss;
				my $end = ($i + 1) * $genebinspan + $tss;
				if ( ($start - int($start)) != 0.0 )
                                {
                                        $start = int ($start + 1);
				}
				if ($i == ($numgenebins - 1)) 
				{	
					if ( int( $end + 0.5 ) != $tts )
					{
						die "for acc: $acc last gene bin boundary is $end, but tts is $tts\n";
					}
					elsif ( $end <= $tts )
					{
						$end += 1; # add one to end to make sure $tts is included in last bin
					}
				}

                                my $plus_cpds = 0;
                                my $plus_dipys = 0;
                                my $minus_cpds = 0;
                                my $minus_dipys = 0;
                                for ( my $j = $start; $j < $end; $j++ )
                                {
                                        if ( exists $plusreads{$j} )
                                        {
                                                $plus_cpds += $plusreads{$j};
                                                $plus_dipys++;
                                        }       
                                        if ( exists $minusreads{$j} )
                                        {
                                                $minus_cpds += $minusreads{$j};
                                                $minus_dipys++;
                                        }
                                }
                                push @cpd_nts, $plus_cpds;
                                push @dipy_nts, $plus_dipys;

                                push @cpd_ts, $minus_cpds;
                                push @dipy_ts, $minus_dipys;
			}

                        # calculate read sums (NMP and Purine bkgd) for terminator
			for ( my $i = 0; $i < $numflankbins; $i++)
                        {
                                my $start = $tts + ( $i * $interval) + 1; # Start one base beyond TTS
                                my $end = $tts + ($i + 1) * $interval;
                                my $plus_cpds = 0;
                                my $plus_dipys = 0;
                                my $minus_cpds = 0;
                                my $minus_dipys = 0;
                                for ( my $j = $start; $j <= $end; $j++ )	# <= end because includes end nucleotide
                                {
                                        if ( exists $plusreads{$j} )
                                        {
                                                $plus_cpds += $plusreads{$j};
                                                $plus_dipys++;
                                        }
                                        if ( exists $minusreads{$j} )
                                        {
                                                $minus_cpds += $minusreads{$j};
                                                $minus_dipys++;
                                        }
                                }
                                push @cpd_nts, $plus_cpds;
                                push @dipy_nts, $plus_dipys;

                                push @cpd_ts, $minus_cpds;
                                push @dipy_ts, $minus_dipys;
                        }

		}
		elsif ( $strand{$acc} eq "-" )
		{
			if ( $tss < $tts )
			{
				die "Gene $acc is listed on minus strand, but the $tss is less than $tts\n";
			}
			# calculate read sums (NMP and Purine bkgd) for promoter
                        for ( my $i = 0; $i < $numflankbins; $i++)
                        {
                                my $end = $tss + $flank_offset - ( $i * $interval); 
                                my $start = $tss + $flank_offset - (($i + 1) * $interval) + 1; # add 1 because 
                                my $plus_cpds = 0;
                                my $plus_dipys = 0;
                                my $minus_cpds = 0;
                                my $minus_dipys = 0;
                                for ( my $j = $start; $j <= $end; $j++ ) # <= end because includes end nucleotide
                                {
                                        if ( exists $plusreads{$j} )
                                        {
                                                $plus_cpds += $plusreads{$j};
                                                $plus_dipys++;
                                        }
                                        if ( exists $minusreads{$j} )
                                        {
                                                $minus_cpds += $minusreads{$j};
                                                $minus_dipys++;
                                        }
                                }

				push @cpd_ts, $plus_cpds;
				push @dipy_ts, $plus_dipys;

				push @cpd_nts, $minus_cpds;
				push @dipy_nts, $minus_dipys;
			}
				
			# calculate read sums for coding region
			for (my $i = 0; $i < $numgenebins; $i++ )
			{
				my $end = $tss - ($i * $genebinspan);
				my $start = $tss - (($i + 1) * $genebinspan);

                                if ($i == ($numgenebins - 1))
                                {
                                        if ( int( $start + 0.5 ) != $tts )
                                        {
                                                die "for acc: $acc last gene bin boundary is $start, but tts is $tts\n";
                                        }
                                        else
                                        {
                                               	$start = $tts;  # start last bin at tts to include tts in coding region
                                        }
                                }
				else
				{
					$start = int ($start + 1)
				}

                                my $plus_cpds = 0;
                                my $plus_dipys = 0;
                                my $minus_cpds = 0;
                                my $minus_dipys = 0;
                                for ( my $j = $start; $j <= $end; $j++ ) # $end nuc is included in the count (e.g., tss)
                                {
                                        if ( exists $plusreads{$j} )
                                        {
                                                $plus_cpds += $plusreads{$j};
                                                $plus_dipys++;
                                        }       
                                        if ( exists $minusreads{$j} )
                                        {
                                                $minus_cpds += $minusreads{$j};
                                                $minus_dipys++;
                                        }
                                }
                                push @cpd_ts, $plus_cpds;
                                push @dipy_ts, $plus_dipys;

                                push @cpd_nts, $minus_cpds;
                                push @dipy_nts, $minus_dipys;
			}

                        # calculate read sums (NMP and Purine bkgd) for terminator
                        for ( my $i = 0; $i < $numflankbins; $i++)
                        {
                                my $end = $tts - ( $i * $interval); 
                                my $start = $tts - ($i + 1) * $interval;
                                my $plus_cpds = 0;
                                my $plus_dipys = 0;
                                my $minus_cpds = 0;
                                my $minus_dipys = 0;
                                for ( my $j = $start; $j < $end; $j++ ) # < end because does not includes end nucleotide (e.g., tts)
                                {
                                        if ( exists $plusreads{$j} )
                                        {
                                                $plus_cpds += $plusreads{$j};
                                                $plus_dipys++;
                                        }
                                        if ( exists $minusreads{$j} )
                                        {
                                                $minus_cpds += $minusreads{$j};
                                                $minus_dipys++;
                                        }
                                }
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
