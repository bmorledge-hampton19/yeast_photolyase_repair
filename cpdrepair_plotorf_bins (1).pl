#!/usr/bin/perl

use strict;
use warnings;

use lib '../';
use GeneCoord;
use CPDReadValues;

# ask for probe filename to analyze
print STDERR "Enter filename of plus strand reads for repair (2hr) timepoint\n";
my $repairplusfile = <STDIN>;
chomp($repairplusfile);

print STDERR "Enter filename of minus strand reads for repair (2hr) timepoint\n";
my $repairminusfile = <STDIN>;
chomp($repairminusfile);

print STDERR "Enter filename of plus strand reads for 0hr timepoint\n";
my $plusfile = <STDIN>;
chomp($plusfile);

print STDERR "Enter filename of minus strand reads for 0hr timepoint\n";
my $minusfile = <STDIN>;
chomp($minusfile);

print STDERR "Loading Gene coordinates\n";
my $genes = GeneCoord->new();
print STDERR "Loading 0hr files\n";
my $zeroreads = CPDReadValues->new($plusfile, $minusfile);

print STDERR "Loading repair (2hr) files\n";
my $repaireads = CPDReadValues->new($repairplusfile, $repairminusfile);

#number of bins to compute
my $numgenebins = 6;
my $numflankbins = 3;

my %chromosomes = $genes->get_chromosomes();
my %trxstart = $genes->get_tss();
my %trxend = $genes->get_tts();
my %strand = $genes->get_strand();

my $interval = 167; # size of flanking interval bins
my $flank_offset = $numflankbins * $interval;

#print header
print "Fraction Lesions Remaining\tNum gene bins=$numgenebins\tData from files - 0hr: $plusfile\t$minusfile\t - repair: $repairplusfile\t$repairminusfile\n";
print "YORF";

for (my $i = 0; $i < $numflankbins; $i++ )
{
        my $start = 0 - $flank_offset + ( $i * $interval);
        my $end = 0 - $flank_offset + (($i + 1) * $interval) - 1;
        print "\t(TS) Promoter ($start to $end)";
}
for ( my $i = 1; $i <= $numgenebins; $i++ )
{
        print "\t(TS) Coding bin $i";
}
for ( my $i = 0; $i < $numflankbins; $i++ )
{
        my $start = ( $i * $interval) + 1;
        my $end = ($i + 1) * $interval;
        print "\t(TS) Terminator ($start to $end)";
}
for (my $i = 0; $i < $numflankbins; $i++ )
{
        my $start = 0 - $flank_offset + ( $i * $interval);
        my $end = 0 - $flank_offset + (($i + 1) * $interval) - 1;
        print "\t(NTS) Promoter ($start to $end)";
}
for ( my $i = 1; $i <= $numgenebins; $i++ )
{
        print "\t(NTS) Coding bin $i";
}
for ( my $i = 0; $i < $numflankbins; $i++ )
{
        my $start = ( $i * $interval) + 1;
        my $end = ($i + 1) * $interval;
        print "\t(NTS) Terminator ($start to $end)";
}
print "\n";

foreach my $chr (sort keys %chromosomes)
{
	print STDERR "Starting $chr\n";
	my %zeroplusreads = $zeroreads->get_plus_reads_for_chromosome($chr);
	my $num_zeroplusreads = scalar keys %zeroplusreads;
	my %zerominusreads = $zeroreads->get_minus_reads_for_chromosome($chr);
	my $num_zerominusreads = scalar keys %zerominusreads;

        my %repairplusreads = $repaireads->get_plus_reads_for_chromosome($chr);
        my $num_repairplusreads = scalar keys %repairplusreads;
        my %repairminusreads = $repaireads->get_minus_reads_for_chromosome($chr);
        my $num_repairminusreads = scalar keys %repairminusreads;
	print STDERR "$chr 0hr reads: $num_zeroplusreads plus reads and $num_zerominusreads minus reads; repair: $num_repairplusreads plus reads and $num_repairminusreads minus reads\n";
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
		
		my @fracremain_cpd_ts = ();
		my @fracremain_cpd_nts = ();
		# calculate read sums (CPD and DIPY bkgd) for gene

		if ( $strand{$acc} eq "+" )
		{
			for ( my $i = 0; $i < $numflankbins; $i++)
			{
				my $start = $tss - $flank_offset + ( $i * $interval);
				my $end = $tss - $flank_offset + (($i + 1) * $interval) - 1; # -1 because we don't include end nuc
 
                        	my $zero_plus_cpds = 0;
                        	my $repair_plus_cpds = 0;
                        	my $zero_minus_cpds = 0;
                        	my $repair_minus_cpds = 0;
				for ( my $j = $start; $j <= $end; $j++ )
				{
	                                if ( exists $zeroplusreads{$j} )
        	                        {
                	                        $zero_plus_cpds += $zeroplusreads{$j};
                        	        }
                                	if ( exists $zerominusreads{$j} )
	                                {
	                                        $zero_minus_cpds += $zerominusreads{$j};
	                                }
	
	                                if ( exists $repairplusreads{$j} )
	                                {
	                                        $repair_plus_cpds += $repairplusreads{$j};
	                                }
	                                if ( exists $repairminusreads{$j} )
	                                {
	                                        $repair_minus_cpds += $repairminusreads{$j};
	                                }
	
				}
	                        my $plus_remaining = "";
	                        if ( $zero_plus_cpds > 0 )
	                        {
	                                $plus_remaining = 1.0 * $repair_plus_cpds / $zero_plus_cpds;
	                        }
	                        my $minus_remaining = "";
	                        if ( $zero_minus_cpds > 0 )
	                        {
	                                $minus_remaining = 1.0 * $repair_minus_cpds / $zero_minus_cpds;
	                        }
	
	                        push @fracremain_cpd_nts, $plus_remaining;
	
                                push @fracremain_cpd_ts, $minus_remaining;
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

                                my $zero_plus_cpds = 0;
                                my $repair_plus_cpds = 0;
                                my $zero_minus_cpds = 0;
                                my $repair_minus_cpds = 0;
                                for ( my $j = $start; $j < $end; $j++ )
                                {
                                        if ( exists $zeroplusreads{$j} )
                                        {
                                                $zero_plus_cpds += $zeroplusreads{$j};
                                        }
                                        if ( exists $zerominusreads{$j} )
                                        {
                                                $zero_minus_cpds += $zerominusreads{$j};
                                        }

                                        if ( exists $repairplusreads{$j} )
                                        {
                                                $repair_plus_cpds += $repairplusreads{$j};
                                        }
                                        if ( exists $repairminusreads{$j} )
                                        {
                                                $repair_minus_cpds += $repairminusreads{$j};
                                        }

                                }       
                                my $plus_remaining = "";
                                if ( $zero_plus_cpds > 0 )
                                {
                                        $plus_remaining = 1.0 * $repair_plus_cpds / $zero_plus_cpds;
                                }
                                my $minus_remaining = "";
                                if ( $zero_minus_cpds > 0 )
                                {
                                        $minus_remaining = 1.0 * $repair_minus_cpds / $zero_minus_cpds;
                                }
                                        
                                push @fracremain_cpd_nts, $plus_remaining;

                                push @fracremain_cpd_ts, $minus_remaining;
			}

			for ( my $i = 0; $i < $numflankbins; $i++)
                        {
                                my $start = $tts + ( $i * $interval) + 1; # Start one base beyond TTS
                                my $end = $tts + ($i + 1) * $interval;
                                my $zero_plus_cpds = 0;
                                my $repair_plus_cpds = 0;
                                my $zero_minus_cpds = 0;
                                my $repair_minus_cpds = 0;
                                for ( my $j = $start; $j <= $end; $j++ )
                                {
                                        if ( exists $zeroplusreads{$j} )
                                        {
                                                $zero_plus_cpds += $zeroplusreads{$j};
                                        }
                                        if ( exists $zerominusreads{$j} )
                                        {
                                                $zero_minus_cpds += $zerominusreads{$j};
                                        }

                                        if ( exists $repairplusreads{$j} )
                                        {
                                                $repair_plus_cpds += $repairplusreads{$j};
                                        }
                                        if ( exists $repairminusreads{$j} )
                                        {
                                                $repair_minus_cpds += $repairminusreads{$j};
                                        }

                                }       
                                my $plus_remaining = "";
                                if ( $zero_plus_cpds > 0 )
                                {
                                        $plus_remaining = 1.0 * $repair_plus_cpds / $zero_plus_cpds;
                                }
                                my $minus_remaining = "";
                                if ( $zero_minus_cpds > 0 )
                                {
                                        $minus_remaining = 1.0 * $repair_minus_cpds / $zero_minus_cpds;
                                }
                                        
                                push @fracremain_cpd_nts, $plus_remaining;

                                push @fracremain_cpd_ts, $minus_remaining;
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
                                my $zero_plus_cpds = 0;
                                my $repair_plus_cpds = 0;
                                my $zero_minus_cpds = 0;
                                my $repair_minus_cpds = 0;
                                for ( my $j = $start; $j <= $end; $j++ )
                                {
                                        if ( exists $zeroplusreads{$j} )
                                        {
                                                $zero_plus_cpds += $zeroplusreads{$j};
                                        }
                                        if ( exists $zerominusreads{$j} )
                                        {
                                                $zero_minus_cpds += $zerominusreads{$j};
                                        }

                                        if ( exists $repairplusreads{$j} )
                                        {
                                                $repair_plus_cpds += $repairplusreads{$j};
                                        }
                                        if ( exists $repairminusreads{$j} )
                                        {
                                                $repair_minus_cpds += $repairminusreads{$j};
                                        }

                                }       
                                my $plus_remaining = "";
                                if ( $zero_plus_cpds > 0 )
                                {
                                        $plus_remaining = 1.0 * $repair_plus_cpds / $zero_plus_cpds;
                                }
                                my $minus_remaining = "";
                                if ( $zero_minus_cpds > 0 )
                                {
                                        $minus_remaining = 1.0 * $repair_minus_cpds / $zero_minus_cpds;
                                }
                                        
                                push @fracremain_cpd_ts, $plus_remaining;

                                push @fracremain_cpd_nts, $minus_remaining;
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
                                my $zero_plus_cpds = 0;
                                my $repair_plus_cpds = 0;
                                my $zero_minus_cpds = 0;
                                my $repair_minus_cpds = 0;
                                for ( my $j = $start; $j <= $end; $j++ )
                                {
                                        if ( exists $zeroplusreads{$j} )
                                        {
                                                $zero_plus_cpds += $zeroplusreads{$j};
                                        }
                                        if ( exists $zerominusreads{$j} )
                                        {
                                                $zero_minus_cpds += $zerominusreads{$j};
                                        }

                                        if ( exists $repairplusreads{$j} )
                                        {
                                                $repair_plus_cpds += $repairplusreads{$j};
                                        }
                                        if ( exists $repairminusreads{$j} )
                                        {
                                                $repair_minus_cpds += $repairminusreads{$j};
                                        }

                                }
                                my $plus_remaining = "";
                                if ( $zero_plus_cpds > 0 )
                                {
                                        $plus_remaining = 1.0 * $repair_plus_cpds / $zero_plus_cpds;
                                }
                                my $minus_remaining = "";
                                if ( $zero_minus_cpds > 0 )
                                {
                                        $minus_remaining = 1.0 * $repair_minus_cpds / $zero_minus_cpds;
                                }

                                push @fracremain_cpd_ts, $plus_remaining;

                                push @fracremain_cpd_nts, $minus_remaining;
			}

                        # calculate read sums (NMP and Purine bkgd) for terminator
                        for ( my $i = 0; $i < $numflankbins; $i++)
                        {
                                my $end = $tts - ( $i * $interval); 
                                my $start = $tts - ($i + 1) * $interval;
                                my $zero_plus_cpds = 0;
                                my $repair_plus_cpds = 0;
                                my $zero_minus_cpds = 0;
                                my $repair_minus_cpds = 0;
                                for ( my $j = $start; $j < $end; $j++ )
                                {
                                        if ( exists $zeroplusreads{$j} )
                                        {
                                                $zero_plus_cpds += $zeroplusreads{$j};
                                        }
                                        if ( exists $zerominusreads{$j} )
                                        {
                                                $zero_minus_cpds += $zerominusreads{$j};
                                        }

                                        if ( exists $repairplusreads{$j} )
                                        {
                                                $repair_plus_cpds += $repairplusreads{$j};
                                        }
                                        if ( exists $repairminusreads{$j} )
                                        {
                                                $repair_minus_cpds += $repairminusreads{$j};
                                        }

                                }
                                my $plus_remaining = "";
                                if ( $zero_plus_cpds > 0 )
                                {
                                        $plus_remaining = 1.0 * $repair_plus_cpds / $zero_plus_cpds;
                                }
                                my $minus_remaining = "";
                                if ( $zero_minus_cpds > 0 )
                                {
                                        $minus_remaining = 1.0 * $repair_minus_cpds / $zero_minus_cpds;
                                }

                                push @fracremain_cpd_ts, $plus_remaining;

                                push @fracremain_cpd_nts, $minus_remaining;
                        }

		}
                if ( scalar @fracremain_cpd_ts != scalar @fracremain_cpd_nts )
                {
                        die "Arrays are of different sizes!\n";
                }

                # bin values

                # print Fraction of CPD's remaining ( unrepaired ) for acc
                print "$acc";

                foreach my $val (@fracremain_cpd_ts)
                {
                        print "\t$val";
                }
                foreach my $val (@fracremain_cpd_nts)
                {
                        print "\t$val";
                }
                print "\n";
        }

}

