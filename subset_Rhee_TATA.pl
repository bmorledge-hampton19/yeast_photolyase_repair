use strict;
use warnings;

print STDERR "What is name of TFBS file\n";
my $file = <STDIN>;
my $dipyout = $file;
$dipyout =~ s/\.fa/_dipy\.bed/ || die "wrong file format\n";
my $notdipy = $file;
open ( DIPY, ">$dipyout" ) || die "couldn't open file\n";
$notdipy =~ s/\.fa/_NOTdipy\.bed/ || die "wrong file format\n";
open ( NOT, ">$notdipy" ) || die "couldn't open file\n";
open ( INPUT, "$file" ) || die "couldn't open file\n";

while ( <INPUT> )
{
	chomp $_;
	my $format = "";
	if ( $_ =~ /^>(chr[XIV]+):([0-9]+)-([0-9]+)\(([+\-])/ )
	{
		$format = "$1\t$2\t$3\t$4\n";
	}
	else
	{	die "Weird line: $_\n";	}
	my $line = <INPUT>;
	chomp $line;
	my $sub = substr $line, 3, 3;
	if ( $sub eq "AAA" )
	{
		print DIPY $format;
	}
	else
	{
		print NOT $format;
	}

}
