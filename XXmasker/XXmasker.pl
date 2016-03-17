#!/usr/bin/env perl

=head1 NAME

XXmakser.pl - mask homologous regions, repeats and low complexity regions

=head1 SYNOPSIS

XXmasker.pl [options] file

=head1 DESCRIPTION

XXmasker

=head1 OPTIONS

=over 8

=item B<--blastdir DIR>

Directory containing the BLAST executables. Default:
/cluster/toolkit/production/blast/bin, then try $PATH, then fail.

=item B<-c N, --cpu N>

Number of threads for BLAST. Default: one thread per available CPU.

=item B<-d FILE, --database FILE>

Database for profile construction.
Default: /cluster/databases/standard/nre70

=item B<--debug>

Show debug information.

=item B<-h, -?, --help>

Print this message and exit.

=item B<--man>

Show full man page.

=item B<-o FILE, --outfile FILE>

Output file. Default: same name as input file with suffix '_hf' after
the basename.

=item B<--cutoff N>

cutoff for homology masking. Default 1e-10.

=back

=head1 DESCRIPTION

Compare a set of sequences to mask homologous regions by
iteratively adding them to a blast database.
First, do one iteration of PSI-BLAST to construct a profile.
Then, run one iteration of BLAST against the database of already
masked sequences. Hit regions are N'ed out.

=cut

use warnings;
use strict;

use lib "lib";

use MiscHelpers;

use Bio::SeqIO;
use Bio::AlignIO;
use File::Basename;
use File::Temp;
use File::Spec;
use Getopt::Long;
use POSIX qw/floor/;

# set global verbosity
Bio::Root::Root::verbose(-1);

sub pushAlignment($$;);
sub lowComplexityFilter($$;);
sub repeatFilter($$;);
sub calcN_frac($;);
sub formatSequence($;);
sub maskShortSegments($$;);

### normally no changes below

# process command line options
my $blastdir      = "";
my $cpu           = 0;
my $debug         = 0;
my $help          = 0;
my $jobspath      = "";
my $man           = 0;
my $outfile       = "";
my $retainedfile  = "";
my $cutoff        = 1e-10;
my $multAli       = 0;
my $lowComplexity = 1;
my $repeatFilter  = 1;
my $backwards     = 0;
my $batch 		  = 0;

print STDERR "====================\n";
print STDERR "==  run XXmasker  ==\n";
print STDERR "====================\n";

if( !GetOptions(
	'blastdir=s' => \$blastdir,
	'c|cpu=i'    => \$cpu,
	'debug!'     => \$debug,
	'multAli'    => \$multAli,
	'h|help|?'   => \$help,
	'jobspath=s' => \$jobspath,
	'man'        => \$man,
	'outfile=s'  => \$outfile,
	'cutoff=f'   => \$cutoff,
	'batch'		 => \$batch
)){
	print STDERR "invalid parameters inserted!!!";
	exit(-1);
}

#if ( @ARGV > 1 ) { pod2usage( -verbose => 1, -exitval => 2 ); }
#if ($help)       { pod2usage( -verbose => 1, -exitval => 0 ); }
#if ($man)        { pod2usage( -verbose => 2, -exitval => 0 ); }

unless ($blastdir) {
	my $blast = `which blastall`;
	#print "blastall: $blast\n";
	if ($blast) {
		( my $basename, $blastdir, my $extension ) =
		  fileparse( $blast, qr/\.[^.]*?/ );
	}else {
		#pod2usage(
		#	-message => "Unable to locate BLAST executables",
		#	-verbose => 1,
		#	-exitval => 1
		#);
		print STDERR "\nUnable to locate BLAST executables. Install blast2!!! \n\n";
		exit(-1);
	}
}

$cpu = $cpu ? $cpu : getNumberOfCPUs();

my $infile = $ARGV[0];
if ( !$infile ) { print STDERR "\nno input file given\n\n"; exit(-1); }
if ( !-e $infile ) { print STDERR "\nUnable to open $infile\n\n"; exit(-1);}
if ( !-r $infile ) { print STDERR "\n$infile is not readable\n\n"; exit(-1);}

unless ($outfile) {
	( my $basename, my $path, my $extension ) =
	  fileparse( $infile, qr/\.[^.]*?/ );
	$outfile = File::Spec->catfile( $path, $basename . "_hf" . $extension );
}

if ($debug) {
	print "Using options:\n";
	print "\tinfile = $infile\n";
	print "\tblastdir = $blastdir\n";
	print "\tcpu = $cpu\n";
	print "\toutfile = $outfile\n";
	print "\tmultAli = $multAli\n";
}

if ($debug) {
	print "Command line parameters after option processing:\n";
	for ( my $i = 0 ; $i < @ARGV ; ++$i ) {
		print "\t$i: $ARGV[$i]\n";
	}
}

my %tmpopts;
#if ($debug) { $tmpopts{CLEANUP} = 0; }

( my $inbasename, my $path, my $ext ) = fileparse( $infile, qr/\.[^.]*?/ );

my $tmpdir = File::Temp->newdir("/tmp/hf_${inbasename}_XXXXX", %tmpopts);
#print "tmpDir: $tmpdir\n";

###
### split input sequences in multiple files
###
open(INPUT, "$infile" ) || die "couldn't open the file! $infile";
my @input = <INPUT>;
close(INPUT);

if(!$multAli){ 
	for(my $i = 0; $i < @input; $i++){
		if($input[$i] =~ /\/\//){
			print STDERR "\nError: input file seems to be in multiple alignment format !\n";
			print STDERR "\tuse option --format MFASTA\n\n";
			exit(-1);
		}
	}
}

my @alignments;
for(my $i = 0; $i < @input; $i++){
	my $line = $input[$i];
	chomp $line;
	if ( $line =~ />(.*)/ ) {
		my $TF = ( split( /\s+/, $1 ) )[0];
		$line =~ s/\s+/;/g;

		#my $inputFileName = "$tmpdir/$TF.fas";
		my $inputFileName = "$tmpdir/tmp.fas";
		open( OUT, ">$inputFileName" ) || die "couldn't open file $inputFileName for writing";
		print OUT $line . "\n";		
		#print $line . "\n";
		my $seq = "";
		$i++;
		while($i < @input && $input[$i] !~ />/){
			$line = $input[$i];
			chomp $line;
			$seq .= $line;			
			$i++;
		}
		$i--;		
		#print "$seq\n";
		my $N_frac = calcN_frac( $seq );				

		print OUT $seq;
		if ($multAli) {
			$i++;
			while($i<@input && $input[$i] !~ /\/\// ){
				$line = $input[$i];
				chomp $line;
				if($line =~ />/){ 
					print OUT "\n" . $line . "\n"; 
				}else{
					print OUT $line;
				}
				$i++;
			}						
		}
		print OUT "\n";
		close(OUT);			

		## push alignment from inputFile in alignments array
		pushAlignment( $inputFileName, \@alignments ) if($N_frac < 1);
		#pushAlignment( $inputFileName, \@alignments );
		unlink $inputFileName;
	}
}

if ($debug) {
	print "\nRead "
	  . scalar(@alignments)
	  . " alignments"
	  . ( @alignments == 1 ? "" : "s" )
	  . " from $infile:\n";

# 	foreach my $aln (@alignments) {
# 		my $s = ($$aln->each_seq)[0];
# 		print ">". $s->id . "\n" . $s->seq . "\n";
# 		
# 		$s = ($$aln->each_seq)[1];
# 		print ">". $s->id . "\n" . $s->seq . "\n";
# 		
# 		print "\/\/\n";
# 	}
# 	exit(-1); 	
}


# keep track how many residues are X'ed out for statistics
my $resTotal  = 0;    # number of residues in all sequences
my $resMasked = 0;

# first, do a dummy search against a database with all sequences
# to determine effective database length
open DB, ">$tmpdir/dummy_db.fasta" or die "Can't create files in $tmpdir\n";
foreach my $aln (@alignments) {
	my $s = ( $$aln->each_seq )[0];
	print DB ">" . $s->id . "\n" . $s->seq . "\n";
}
close DB;

open QU, ">$tmpdir/dummy_query.fasta" or die "Can't create files in $tmpdir\n";
my $s = ( ${ $alignments[0] }->each_seq )[0];
print QU ">" . $s->id . "\n" . $s->seq . "\n";
#print ">" . $s->id . "\n" . $s->seq . "\n";
close QU;
#print "$blastdir/formatdb -p F -i $tmpdir/dummy_db.fasta -n $tmpdir/dummy_db\n";
`$blastdir/formatdb -p F -i $tmpdir/dummy_db.fasta -n $tmpdir/dummy_db -l $tmpdir/formatdb.log`;

my $effDBSize;
open BLA,
"$blastdir/blastall -F F -p blastn -a $cpu -i $tmpdir/dummy_query.fasta -d $tmpdir/dummy_db |";
while ( my $line = <BLA> ) {
	next unless $line =~ s/^Effective length of database: //;
	$line =~ s/,//g;
	chomp($line);
	$effDBSize = $line;
	last;
}
close BLA;

if ($debug) {
	print "\nEffective database size set to $effDBSize\n";
}

if ($debug) {
	print "Profile construction completed\n";
	print "Starting DB construction\n";
}

# fasta file basename with sequences for blast database
my $outDB       = "$tmpdir/filtered_db";
my $outDB_fasta = "$tmpdir/filtered_db.fasta";
open OUTDB, ">$outDB_fasta" or die "$!\n";

###
# main loop
# For each sequence, run blast against database so far, mask hits and
# add masked sequence to database.
###

my $allOutSeqs = "";

for ( my $i = 0 ; $i < @alignments ; ++$i ) {
	## print Progress Bar
	if ( !$batch && !$debug && $i % 100 ) {
		my $width   = 31;
		my $percent = $i / $#alignments;
		print "\r mask sequences:            [";
		print "=" for ( 1 .. floor( $percent * $width ) );
		print " " for ( 1 .. ( $width - floor( $percent * $width ) ) );
		printf "] %.2f%%", $percent * 100;
	}
	my $s = ( ${ $alignments[$i] }->each_seq )[0];
	# write file with only first sequence of alignment
	open FAS, ">$tmpdir/$i.fasta";
	print FAS ">" . $s->id . "\n" . $s->seq . "\n";
	close FAS;

	$resTotal += length( $s->seq );

	my $seqMasked = $s->seq;
	$seqMasked =~ s/N/%/g;    # discriminate genuine X's from masking

	if ( -z OUTDB ) {         # first sequence: add to database
		if ($debug) {
			print "Sequence $i (" . $s->id . ") accepted without masking (first sequence)\n";
		}
	} else {
		# now run blast against database so far
		my $cmd = "$blastdir/blastall -F \"m D\" -p blastn -a $cpu -e $cutoff -z $effDBSize -i $tmpdir/$i.fasta -d $outDB > $tmpdir/$i.bla2";
		#my $cmd = "$blastdir/blastall -F F -p blastn -a $cpu -e $cutoff -z $effDBSize -i $tmpdir/$i.fasta -d $outDB > $tmpdir/$i.bla2";
		sysExec($cmd);

		# collect hit regions from blast output
		open BLA, "$tmpdir/$i.bla2" or die "Can't open $i.bla2\n";
		my @hits;
		my @hitnames;
		my $sbjctName;
		my $firstQueryLine;
		my $lastQueryLine;
		while ( my $line = <BLA> ) {
			if ( $line =~ /^ Score/ || eof(BLA) ) {

				# if there was a hit so far, push it to hit list
				if ($firstQueryLine) {
					my @f = split( '\s+', $firstQueryLine );
					my @l = split( '\s+', $lastQueryLine );
					push( @hits,     "$f[1]-$l[3]" );
					push( @hitnames, "$f[1]-$l[3] ($sbjctName)" );
					$firstQueryLine = "";
				}
			}
			elsif ( $line =~ /^Query:/ ) {
				$firstQueryLine = $firstQueryLine ? $firstQueryLine : $line;
				$lastQueryLine = $line;
			}
			if ( $line =~ /^>/ ) {
				$sbjctName = substr( $line, 1 );
				chomp($sbjctName);
			}
		}
		if ( $debug && @hits ) {
			print "\n" . "-" x 80 . "\n";
			print "Running BLASTN against masked DB for sequence $i (" . $s->id . ")\n$cmd\n";
			print scalar(@hits) . " hit" . ( @hits > 1 ? "s" : "" ) . ": " . join( ', ', @hitnames ) . "\n";
		}
		else {
			#print "No hits\n";
		}
		close BLA;
		unlink "$tmpdir/$i.fasta";
		unlink "$tmpdir/$i.bla2";
		
		$seqMasked = lowComplexityFilter( $i, $seqMasked ) if ($lowComplexity);		
		$seqMasked = repeatFilter( $i, $seqMasked ) if ($repeatFilter);
		
		for my $hit (@hits) {
			( my $start, my $end ) = split( '-', $hit );
			--$start;
			--$end;
			my $len = $end - $start + 1;	
					
			substr( $seqMasked, $start, $len ) = 'N' x $len;
		}
		
		#$seqMasked = maskShortSegments( $min_length, $seqMasked );
		
		# do statistics
		my $numResMasked = 0;
		$numResMasked++ while $seqMasked =~ /N/g;
		$resMasked += $numResMasked;
		
		if ( $debug && @hits ) {
			print "  masked sequence (" . length($seqMasked) . " residues, $numResMasked masked by homology):\n";
			print formatSequence($seqMasked) . "\n";
		}
	}

	$seqMasked =~ s/%/N/g;    # restore original N's

	my $seqFormatted = formatSequence($seqMasked);
	$seqFormatted =~ s/ //g;
	$seqFormatted =~ s/\n//g;

	my $idLine = $s->id;
	$idLine =~ s/;/ /g;	
	#my $N_frac = calcN_frac($seqFormatted);
   #if($N_frac < 0.95){
		#my @idsplit = split(/;/, $s->id);
		my $outSeq = ">$idLine\n$seqFormatted\n";

		#print $s->id . " " . $idsplit[0] . " " . $idsplit[1] . "\n";
		#exit(-1);
		#my $outSeq = ">".$s->id ."\n$seqFormatted\n";

		# write unmasked sequences of alignment into results file
		for ( my $j = 1 ; $j < ${ $alignments[$i] }->each_seq ; $j++ ) {
			my $s_ali = ( ${ $alignments[$i] }->each_seq )[$j];
			$outSeq .= ">" . $s_ali->id . "\n" . $s_ali->seq . "\n";
		}
		if ($multAli) {
			$outSeq .= "//\n";
		}
		if($backwards == 1){
			$allOutSeqs = $outSeq . $allOutSeqs;
		}else{
			$allOutSeqs = $allOutSeqs . $outSeq;
		}
	#}
	# add original sequence to database and format it for next round
	print OUTDB ">" . $s->id . "\n" . $s->seq . "\n";

	#print "formatting database $outDB\n" if $debug;
	system("$blastdir/formatdb -p F -i $outDB_fasta -n $outDB -l $tmpdir/formatdb.log") == 0 || die "error formatting blast db\n";

}

close OUTDB;

print "\n" if(!$batch);
print "XXmasker: $resMasked of $resTotal residues masked because of homology\n\n";
#print $allOutSeqs;
open OUT, ">$outfile" or die "$!: $outfile\n";
print OUT $allOutSeqs;
close OUT;
unlink "error.log";

# clean up tmp directory
# system("rm -rf $tmpdir") == 0 or die "$!: $tmpdir\n";

sub pushAlignment($$;) {
	my $inputFileName = $_[0];
	my $alignments    = $_[1];
	
	my $aln = Bio::AlignIO->new( -file => $inputFileName, -format => 'fasta' ) ->next_aln();
	if($backwards == 1){
		unshift @$alignments, \$aln;
	}else{
		push @$alignments, \$aln;
	}
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
	while ( $seq =~ /(^|X)([^X]{1,$thresh})(?:X|$)/ ) {
		substr( $seq, length($`) + length($1), length($2) ) = 'X' x length($2);
	}
	return $seq;
}

sub lowComplexityFilter($$;) {
	my $id           = $_[0];
	my @seq          = split( //, $_[1] );
	my $minRepLength = 50;

	my $seqlen    = $#seq;
	my $repLength = 2;
	my @nucs      = ( $seq[0], $seq[1] );
	for ( my $i = 2 ; $i <= $seqlen + 1 ; $i++ ) {
		my $base;
		if   ( $i == $seqlen + 1 ) { $base = 0; }
		else                       { $base = $seq[$i] }
		if ( $nucs[0] eq $base || $nucs[1] eq $base ) {
			$repLength++;
		}
		elsif ( $nucs[0] eq $nucs[1] ) {
			$repLength++;
			$nucs[1] = $base;
		}
		else {
			if ( $repLength > $minRepLength && $nucs[0] ne 'N' && $nucs[1] ne 'N' )
			{
				printf "mask low complexity of length %d in %s (%s, %s)\n",
				  $repLength, $id, $nucs[0], $nucs[1]  if ($debug);
				for ( my $j = 1 ; $j <= $repLength ; $j++ ) {
					$seq[ $i - $j ] = 'N';
				}
			}
			else {
				$i -= ( $repLength - 1 );
			}
			$nucs[0]   = $seq[ $i++ ];
			$nucs[1]   = $seq[$i];
			$repLength = 2;
		}
	}
	return join( "", @seq );
}

sub repeatFilter($$;) {
	my $id               = $_[0];
	my @seq              = split( //, $_[1] );
	my $minUnitLength    = 3;
	my $maxUnitLength    = 10;
	my $minStretchLength = 50;

	my $seqlen = $#seq;
	for (
		my $unitLength = $minUnitLength ;
		$unitLength <= $maxUnitLength ;
		$unitLength++
	  )
	{
		my @nucs;
		for ( my $j = 0 ; $j < $unitLength ; $j++ ) {
			push( @nucs, $seq[$j] );
		}
		my $repLength = $unitLength;
		for ( my $i = $unitLength + 1 ; $i <= $seqlen + 1 ; $i++ ) {
			my $base;
			if   ( $i == $seqlen + 1 ) { $base = 0; }
			else                       { $base = $seq[$i] }
			if ( $nucs[ $i % $unitLength ] eq $base ) {
				$repLength++;
			}
			else {
				my $N_rep = 0;
				for ( my $j = 0 ; $j < $unitLength ; $j++ ) {
					$N_rep = 1 if ( $nucs[$j] eq 'N' );
				}
				if ( $repLength > $minStretchLength && !$N_rep ) {
					if ($debug) {
						printf "mask repeat of size %d and length %d in %s (",
						  $unitLength, $repLength, $id;
						for ( my $j = 0 ; $j < $unitLength ; $j++ ) {
							print $seq[ $i - $repLength + $j ];
						}
						print ")\n";
					}
					for ( my $j = 1 ; $j <= $repLength ; $j++ ) {
						$seq[ $i - $j ] = 'N';
					}
				}
				else {
					$i -= ( $repLength - 1 );
				}
				for ( my $j = 0 ; $j < $unitLength ; $j++, $i++ ) {
					$nucs[ $i % $unitLength ] = $seq[$i];
				}
				$i--;
				$repLength = $unitLength;
			}
		}
	}
	return join( "", @seq );
}

sub calcN_frac($;){
	return 1 if(length $_[0] == 0);
	my @seq = split( //, $_[0] );
  my $length = $#seq + 1;
  my $N = 0;
	for(my $i=0; $i<@seq; $i++){
	   $N++ if($seq[$i] eq 'N');
	}	
	return $N/$length;
}
