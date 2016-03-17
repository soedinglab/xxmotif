#!/usr/bin/env perl

=head1 NAME

homology_filter.pl - mask homologous regions in a set of protein sequences

=head1 SYNOPSIS

homology_filter.pl [options] file

=head1 OPTIONS

=over 8

=item B<--blastdir DIR>

Directory containing the BLAST executables. Default:
/cluster/toolkit/production/bioprogs/blast/bin, then try $PATH, then fail.

=item B<-c N, --cpu N>

Number of threads for BLAST. Default: 1.

=item B<-d FILE, --database FILE>

FASTA file from which database is to be built.

=item B<--debug>

Show debug information.

=item B<-h, -?, --help>

Print this message and exit.

=item B<--man>

Show full man page.

=item B<--min-length N>

Residue stretcher shorter than or as long as N which lie between
two homology masked regions or between one region and the end of
the sequence will be masked, too, as well as sequences with
total length <=N. Default: 25.

=item B<-o FILE, --outfile FILE>

Output file. Default: same name as input file with suffix '_hf' after
the basename.

=item B<--retain-threshold [0.0,1.0+eps]>

All sequences with a fraction of less than THRESH masked positions
will be included in the filtered set. Default: 0.5.

=back

=head1 DESCRIPTION

Compare a set of ELM sequences to mask homologous regions by
iteratively adding them to a blast database.
First, do one iteration of PSI-BLAST to construct a profile.
Then, run one iteration of BLAST against the database of already
masked sequences. Hit regions are X'ed out.

=cut

use warnings;
use strict;

use Bio::SeqIO;
use Data::Dumper;
use File::Basename;
use File::Copy;
use File::Temp qw(tempdir);
use File::Spec;
use Getopt::Long;
use Pod::Usage;

sub formatSequence($;);
sub maskShortSegments($$;);

### normally no changes below

# process command line options
my $blastdir   = "";
my $cpu        = 1;
my $database   = "";
my $debug      = 0;
my $help       = 0;
my $man        = 0;
my $min_length = 25;
my $outfile    = "";
my $retthresh  = 0.5;
my $dcinfile   = "";    # disorder/conservation, used if available
my $dcoutfile  = "";
my $propinfile   = "";    # NN predictor properties, used if available
my $propoutfile  = "";
my $usedc      = 0;
my $useprop    = 0;

GetOptions(
	'blastdir=s'         => \$blastdir,
	'c|cpu=i'            => \$cpu,
	'd|database=s'       => \$database,
	'debug!'             => \$debug,
	'h|help|?'           => \$help,
	'man'                => \$man,
	'min-length=i'       => \$min_length,
	'outfile=s'          => \$outfile,
	'retain-threshold=f' => \$retthresh,
) or pod2usage(-verbose => 1, -exitval => 2);

if (@ARGV > 1) {
	print STDERR "Too many arguments (" . join('; ', @ARGV) . ")\n";
	pod2usage(-verbose => 1, -exitval => 2);
}

if ($help) {
	pod2usage(-verbose => 1, -exitval => 0);
}

if ($man) {
	pod2usage(-verbose => 2, -exitval => 0);
}

unless ($blastdir) {
	my $default = '/cluster/toolkit/production/bioprogs/blast/bin';
	if (-e $default) {
		$blastdir = $default;
	} else {
		my $bla = `which blastpgp`;
		print "blastpgp: $bla\n" if $debug;
		if ($bla) {
			(my $basename, $blastdir, my $extension) =
			  fileparse($bla, qr/\.[^.]*?/);
		} else {
			pod2usage(
				-message => "Unable to locate BLAST executables",
				-verbose => 1,
				-exitval => 1
			);
		}
	}
}
$database = "/cluster/databases/standard/nre70" unless ($database);

my $infile = $ARGV[0];
if (!$infile) {
	pod2usage(-verbose => 1, exitval => 2);
}
if (!-e $infile) {
	pod2usage("Unable to open $infile");
}
if (!-r $infile) {
	pod2usage("$infile is not readable");
}

{
	(my $basename, my $path, my $extension) = fileparse($infile, qr/\.[^.]*?/);
	unless ($outfile) {
		$outfile = File::Spec->catfile($path, $basename . "_hf" . $extension);
	}
	$dcinfile = File::Spec->catfile($path, $basename . ".dc");
	$propinfile = File::Spec->catfile($path, $basename . ".prop");
	($basename, $path, $extension) = fileparse($outfile, qr/\.[^.]*?/);
	$dcoutfile = File::Spec->catfile($path, $basename . ".dc");
	$propoutfile = File::Spec->catfile($path, $basename . ".prop");
}

if ($debug) {
	print "Homology filter running on: " . `uname -a`. "\n";
	print "Using options:\n";
	print "\tinfile = $infile\n";
	print "\tblastdir = $blastdir\n";
	print "\tcpu = $cpu\n";
	print "\tdatabase = $database\n";
	print "\toutfile = $outfile\n";
	print "\tretainthreshold = $retthresh\n";
	print "\tdcinfile = $dcinfile\n";
	print "\tpropinfile = $propinfile\n";
}

if ($debug) {
	print "Command line parameters after option processing:\n";
	for (my $i = 0; $i < @ARGV; ++$i) {
		print "\t$i: $ARGV[$i]\n";
	}
}

# log file
(my $inbasename, my $path, my $ext) = fileparse($infile, qr/\.[^.]*?/);

###
### Read input sequences
###
my $in = Bio::SeqIO->new(-file => "$infile", '-format' => 'Fasta');
my @sequences;
while (my $seq = $in->next_seq()) {
	my %s = (
		seq => $seq->seq(),
		id  => $seq->id,
	);
	push @sequences, \%s;
}

if ($debug) {
	print "\nRead "
	  . scalar(@sequences)
	  . " sequence"
	  . (@sequences == 1 ? "" : "s")
	  . " from $infile:\n";
	foreach my $s (@sequences) {
		print ">$$s{id}\n$$s{seq}\n";
	}
}

###
### Read disorder/conservation, if available
###
my %disocons;
if (-e $dcinfile) {
	$usedc = 1;
	open DC, $dcinfile;
	my $id;
	while (my $line = <DC>) {
		chomp($line);
		if ($line =~ /^>/) {
			($id = $line) =~ s/^>//;
		} else {
			die "Illegal disorder/conservation file $dcinfile\n";
		}
		my $diso = <DC>;
		chomp($diso);
		$disocons{$id}{diso} = $diso;
		my $cons = <DC>;
		chomp($cons);
		$disocons{$id}{cons} = $cons;
	}
	close DC;
} else {
	warn "$dcinfile does not exist, not generating $dcoutfile\n";
}
my %props;
if (-e $propinfile) {
	$useprop = 1;
	open PROP, $propinfile;
	my $id;
	while (my $line = <PROP>) {
		chomp($line);
		if ($line =~ /^>/) {
			($id = $line) =~ s/^>//;
		} else {
			die "Illegal property file $propinfile\n";
		}
		for (my $i=0; $i<5; ++$i) {
			my $prop = <PROP>;
			chomp($prop);
			$props{$id} .= $prop . "\n";
		}
	}
	close PROP;
} else {
	warn "$propinfile does not exist, not generating $propoutfile\n";
}

###
### change to temporary directory
###
my $tmpdir = tempdir('hf_XXXXX', TMPDIR => 1, CLEANUP => $debug ? 0 : 1);

# format blast database
copy($database, File::Spec->catfile($tmpdir, "in.fasta")) or die $!;
system(
	"'$blastdir/formatdb' -i '$tmpdir/in.fasta' -n '$tmpdir/in_set_db' -l /dev/null"
  ) == 0
  or die "error formatting blast db\n";

# split into single fasta files and prepare blast jobscripts
for (my $i = 0; $i < @sequences; ++$i) {
	open FAS, '>', File::Spec->catfile($tmpdir, $i . ".fasta");
	print FAS ">$sequences[$i]{id}\n$sequences[$i]{seq}\n";
	close FAS;
}

# keep track how many residues are X'ed out for statistics
my $resTotal  = 0;    # number of residues in all sequences
my $resMasked = 0;

###
# main loop
# For each sequence, run blast against database so far, mask hits and
# add masked sequence to database.
###
my $i            = 0;
my $allOutSeqs   = "";
my $allDcOutSeqs = "";
my $allPropOutSeqs = "";
for my $s (@sequences) {
	print "\n" . "-" x 80 . "\n" if $debug;
	$resTotal += length($$s{seq});

	my $seqMasked = $$s{seq};
	$seqMasked =~ s/X/%/g;    # discriminate genuine X's from masking

	my $cmd =
	  "'$blastdir/blastpgp' -a $cpu -e 0.01 -i '$tmpdir/$i.fasta' -d '$tmpdir/in_set_db' >'$tmpdir/$i.bla2' 2>&1";
	print "Running PSI-BLAST against $database for sequence $i ($$s{id}\n$cmd\n"
	  if $debug;
	system($cmd);

	# collect hit regions from blast output
	open BLA, '<', File::Spec->catfile($tmpdir, $i . ".bla2")
	  or die "Can't open '$tmpdir/$i.bla2'\n";
	my @hits;
	my @hitnames;
	my $sbjctName;
	my $firstQueryLine;
	my $lastQueryLine;
	while (my $line = <BLA>) {

		if ($line =~ /^ Score/ || eof(BLA)) {

			# if there was a hit so far, push it to hit list
			if ($firstQueryLine) {
				my @f = split('\s+', $firstQueryLine);
				my @l = split('\s+', $lastQueryLine);
				push(@hits,     "$f[1]-$l[3]");
				push(@hitnames, "$f[1]-$l[3] ($sbjctName)");
				$firstQueryLine = "";
			}
		} elsif ($line =~ /^Query:/) {
			$firstQueryLine = $firstQueryLine ? $firstQueryLine : $line;
			$lastQueryLine = $line;
		}
		if ($line =~ /^>/) {
			$sbjctName = substr($line, 1);
			chomp($sbjctName);
		}
	}
	if ($debug && @hits) {
		print scalar(@hits) . " hit"
		  . (@hits > 1 ? "s" : "") . ": "
		  . join(', ', @hitnames) . "\n";
	} else {
		print "No hits\n" if $debug;
	}
	close BLA;

	if ($debug) {
		print "original sequence (" . length($$s{seq}) . " residues):\n";
		print formatSequence($$s{seq}) . "\n";
	}
	for my $hit (@hits) {
		(my $start, my $end) = split('-', $hit);
		--$start;
		--$end;
		my $len = $end - $start + 1;
		substr($seqMasked, $start, $len) = 'X' x $len;
	}
	$seqMasked = maskShortSegments($min_length, $seqMasked);

	# do statistics
	my $numResMasked = 0;
	$numResMasked++ while $seqMasked =~ /X/g;
	$resMasked += $numResMasked;

	if ($debug && @hits) {
		print "  masked sequence ("
		  . length($seqMasked)
		  . " residues, $numResMasked masked by homology):\n";
		print formatSequence($seqMasked) . "\n";
	}

	# determine fraction of masked residues (by homology only)
	my $maskedcount = 0;
	for (my $i = 0; $i < length($seqMasked); ++$i) {
		++$maskedcount if substr($seqMasked, $i, 1) eq 'X';
	}
	my $maskedFraction = $maskedcount / length($seqMasked);
	$seqMasked =~ s/%/X/g;    # restore original X's

	if ($maskedFraction < $retthresh) {

		# write out masked sequence to fasta file
		my $seqFormatted = formatSequence($seqMasked);
		$seqFormatted =~ s/ //g;
		$allOutSeqs .= ">$$s{id}\n$seqFormatted\n";
		if ($usedc) {
			if (!defined($disocons{$$s{id}})) {
				die "ERROR: no diso/cons for $$s{id}\n";
			}
			$allDcOutSeqs .=
			  ">$$s{id}\n$disocons{$$s{id}}{diso}\n$disocons{$$s{id}}{cons}\n";
		}
		if ($useprop) {
			if (!defined($props{$$s{id}})) {
				die "ERROR: no properties for $$s{id}\n";
			}
			$allPropOutSeqs .= ">$$s{id}\n$props{$$s{id}}";
		}
	} else {
		print "$$s{id} dropped because "
		  . sprintf("%5.1f", $maskedFraction * 100)
		  . "% masked\n"
		  if $debug;
	}

	++$i;
}

print "$resMasked of $resTotal residues masked because of homology\n" if $debug;

open OUT, '>', $outfile or die "`pwd`\nxxx$!: '$outfile'\n";
print OUT $allOutSeqs;
close OUT;

if ($usedc) {
	open OUT, '>', $dcoutfile or die "`pwd`\n$!: '$dcoutfile'\n";
	print OUT $allDcOutSeqs;
	close OUT;
}
if ($useprop) {
	open OUT, '>', $propoutfile or die "`pwd`\n$!: '$propoutfile'\n";
	print OUT $allPropOutSeqs;
	close OUT;
}

sub formatSequence($;) {
	my $seq = shift;
	$seq =~ s/(.{70})/$1\n/g;
	$seq =~ s/(.{10})/$1 /g;
	$seq =~ s/\s+$//;
	return $seq;
}

sub maskShortSegments($$;) {
	my $thresh = $_[0];
	my $seq    = $_[1];
	while ($seq =~ /(^|X)([^X]{1,$thresh})(?:X|$)/) {
		substr($seq, length($`) + length($1), length($2)) = 'X' x length($2);
	}
	return $seq;
}
