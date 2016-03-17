package MiscHelpers;
require Exporter;
use POSIX qw(floor tmpnam);
use List::Util qw(min max);
use Carp;
use Data::Dumper;

our @ISA    = qw(Exporter);
our @EXPORT = (
	'contains',            'getNumberOfCPUs',
	'gnuplotPslatexToEPS', 'histoCount',
	'makeTexString',       'median',
	'movingMedian',        'printHash',
	'sysExec',
);
our @EXPORT_OK = qw();

# returns median element of list
sub median (@) {
	my @list = sort { $a <=> $b } @_;
	my $count = @list;
	my $m;
	if ( $count % 2 ) {
		$m = $list[ floor( $count / 2 ) ];
	} else {
		$m = ( $list[ $count / 2 ] + $list[ $count / 2 - 1 ] ) / 2;
	}
	return $m;
}

# For the given array, compute the moving median with the
# given window size.
# Returns new array with computed values.
sub movingMedian(\@$) {
	if ($_[1]%2==0) {
		warn("Moving median: only odd window sizes make sense. Using window size " . ($_[1]+1) . " instead of $_[1]\n");
	}
	my $WINDOW = int( $_[1] / 2 );
	my @movMed;
	for ( my $i = 0 ; $i < @{ $_[0] } ; ++$i ) {
		my @values;
		for (
			my $j = max( 0, $i - $WINDOW ) ;
			$j <= min( scalar( @{ $_[0] } ) - 1, $i + $WINDOW ) ;
			++$j
		  )
		{
			push( @values, ${ $_[0] }[$j] );
		}
		push( @movMed, median(@values) );
	}
	return @movMed;
}

# increment i-th element of an array, where i is calculated
# by scaling the second parameter with the number of elements in the array
sub histoCount (\@$;) {
	++${ $_[0] }[ min( @{ $_[0] } - 1, floor( @{ $_[0] } * $_[1] ) ) ];
}

sub getNumberOfCPUs () {
	open INFO, "</proc/cpuinfo" or return 0;
	my $num = 0;
	while ( my $line = <INFO> ) {
		if ( $line =~ /^processor/ ) {
			my @cols = split( /\s+/, $line );
			$num = max( $num, $cols[2] );
		}
	}
	++$num;
	return $num;
}

# check if a value is contained in given array
sub contains ($@) {
	my $key = shift;
	for my $s (@_) {
		if ( $key eq $s ) { return 1; }
	}
	return 0;
}

sub printHash($;) {
	my $hash = shift;
	foreach my $k ( keys(%$hash) ) {
		print "$k:\t", $$hash{$k}, "\n";
	}
}

sub sysExec($;) {
	my $cmd = $_[0];
	`$cmd`;
	my $ret = $?;
	if ( $ret != 0 ) {
		my $hostname = `hostname`;
		warn "ERROR executing $cmd (return value: $ret) on $hostname";
	}
	return $ret;
}

#
# create an eps figure from gnuplot input using terminal pslatex
#
sub gnuplotPslatexToEPS($;) {
	my $gpfile = $_[0];
	( my $basename = $gpfile ) =~ s/\.gp$//;
	my $tmpfile = tmpnam();
	( my $tmpdir = $tmpfile ) =~ s/(.*)\/.+/$1/;
	open IN, "gnuplot <$gpfile |"
	  or confess "Can't open gnuplot pipe from input $gpfile\n";
	open TMP, ">$tmpfile" or confess "Can't open temporary file $tmpfile\n";
	print TMP '
\documentclass{article}
\usepackage[a0paper]{geometry}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\pagestyle{empty}
\usepackage{color,sfmath,upgreek,graphicx}
\makeatletter
\fboxsep0mm
\fboxrule0mm
\begin{document}
\small\sffamily\fcolorbox{white}{white}{%
% white box corrects eps bounding box in some strange cases
';

	while ( my $line = <IN> ) {
		$line =~ s/\\endinput//;
		print TMP $line;
	}
	print TMP '}
\end{document}';
	close IN;
	close TMP;
	`cd $tmpdir && latex $tmpfile`;
	`cd $tmpdir && dvips -E -o $basename.eps $tmpfile`;
}

#
# escapes characters which TeX doesn't like in a given string
#
sub makeTexString($;) {
	my $texRegex = $_[0];
	$texRegex =~ s/_/\\_/g;
	$texRegex =~ s/\$/\\\$/g;
	$texRegex =~ s/\{/\\{/g;
	$texRegex =~ s/\}/\\}/g;
	$texRegex =~ s/\^/\\\^{}/g;

	#	$texRegex =~ s/\|/\\\|/g;
	return $texRegex;
}

1;
