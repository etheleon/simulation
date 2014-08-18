#!/usr/bin/env perl

use v5.10;
use strict;

die "usage: $0 chosenGenomes.fna template.fq output abundanceInfo 
		(indel-rate: proportion of indel over all errors; default=0.1)\n" unless $#ARGV>=3;

##################################################
# Part0: Init
say "1. Initializing ...";
##################################################

#Hashes###
my (%seq,%score,%ntrate,%abundance,%range,%globalnt,%gitaxid,%leng);
#	%seq 	taxaID::concatenated sequence of chromosomes
#	%ntrate	nucleotide::the frequencies
#	%abundance taxaID::abundance
#	%range 	taxaID::the lower and upper bounds
#	%globalnt nucleotide::frequency
#	%gitaxid taxid::gi
#	%leng taxid::sequencelength
##########

#Arrays###
my (@taxid,@zipped);						
#	@taxid	stores taxid in decreasing abundances
#	@zipped	stores taxid, lower and upper abundances
##########

#Scalar##
my $totalAbu;						#summed abundance
my $id = 0;						#For adding to fastq header
my $indelrate=0.01;				 	#the INDEl error rate:percentage of errors to be indel
my $outputfile = $ARGV[2]; $outputfile =~ s/.+\///; 	#the outputfile name
##########

#Phred error probabilities
map {$score{$_}=phred($_)} 0..100; 			#Probability 

##################################################
# Part1: read genomes 
say "2. Reading Genomes ...";
##################################################

my $seq = read_fasta_file($ARGV[0]);
%seq =  %$seq;
map { $leng{$_} = length($seq{$_}) } keys %seq;
say "\tStored ",scalar keys %seq, " sequences";

####################################################
## Part2: abundance
say "3. Reading abundances.. ";
###################################################

open my $abundance, '<', $ARGV[3]; #out/sim.0201.out.txt
while(<$abundance>){unless ($. == 1){ 
	my ($abu, $taxid) = (split(/\t/))[1,3];
	$abundance{$taxid} = $abu
}}

##abundance summary
map { $totalAbu += $abundance{$_} } keys %abundance;	#total reads per million abundance 
map { $abundance{$_} = ($abundance{$_}/$totalAbu)} keys %abundance;

#calculate Nucleotide Frequency
map { ntfreq($_) } keys %abundance;

say "\tNucleotide frequencies";
map {say "\t\t$_: $globalnt{$_}"} keys %globalnt;

###taxid - range creation:: TAXID: (lower, upper)
@taxid= sort { $abundance{$b} <=> $abundance{$a} } keys %abundance;	#sort taxid by their abundance
@zipped = map {($taxid[$_], $abundance{$taxid[$_]})} 0 .. $#taxid;

##THINGS TO DO to convert into a mappable function
my $rollingvalue; 
for(my $n=0; $n<=$#zipped; $n+=2)
	{
	    my $taxid = $zipped[$n];	#taxid
		if($n==0){
		    push @{ $range{$taxid} }, 0, $zipped[$n+1];
		    $rollingvalue = $zipped[$n+1];
	    	}else{
		    push @{ $range{$taxid} }, $rollingvalue, $rollingvalue + $zipped[$n+1];
		    $rollingvalue += $zipped[$n+1];
		}
	}

###################################################
## Simulate Shotgun sequence
say "4. Reading quality scores & simulating ...";
###################################################

open my $fastaOutput, 	'>', 	"$ARGV[2]".'.fna' || die $!;				#FASTA OUTPUT
open my $fastqOutputONE, 	'>', 	"$ARGV[2]".'_1'.'.fq'  || die $!;				#FASTQ OUTPUT
open my $fastqOutputTWO, 	'>', 	"$ARGV[2]".'_2'.'.fq'  || die $!;				#FASTQ OUTPUT

#say "reading file: ".$ARGV[1]."_1.filtered.fastq.gz";
my $fastqfileone = $ARGV[1]."_1.filtered.fastq.gz";

#say "reading file: ".$ARGV[1]."_2.filtered.fastq.gz";
my $fastqfiletwo = $ARGV[1]."_2.filtered.fastq.gz";

open my $fastqONE, 	"gzcat $fastqfileone | " || die $!;
open my $fastqTWO, 	"gzcat $fastqfiletwo | " || die $!;

#assumming single line fastq
while (!eof($fastqONE) and !eof($fastqTWO)) {
	<$fastqONE>;<$fastqONE>;<$fastqONE>;	#h1, sequence, h2
	<$fastqTWO>;<$fastqTWO>;<$fastqTWO>;	#h1, sequence, h2

#Process Quality

	my $qualONE 	=  <$fastqONE>;
	chomp $qualONE;
#	say $qualONE;

	my $qualTWO =  <$fastqTWO>;
	chomp $qualTWO;

	my $readlength = length($qualONE);	#template length	#assuming pair1 & 2 have the same read length
	my $qualitystringONE = processQual($qualONE);
	my $qualitystringTWO = processQual($qualTWO);
	
##Choosing taxa and loc to pluck sequence out from

	my $taxaofchoice = choosetaxa(@taxid);
	processPairReads($taxaofchoice, $readlength, $qualitystringONE, $qualitystringTWO, $qualONE,$qualTWO);
}

say "All done. Please check $ARGV[2] for results";
####################################################################################################
#Functions
####################################################################################################

#Calculates error-probability of mutation
sub phred
{
my ($score) = @_;
	return 10**($score/(-10));
}

sub ntfreq
{
my ($taxid) = @_;
    	my %ntrate;
	$ntrate{'a'} += ($seq{$taxid}=~tr/aA/AA/);	#count the number of ATCGs
	$ntrate{'t'} += ($seq{$taxid}=~tr/tT/TT/);
	$ntrate{'g'} += ($seq{$taxid}=~tr/gG/GG/);
	$ntrate{'c'} += ($seq{$taxid}=~tr/cC/CC/);
	map {	$globalnt{$_} += ($ntrate{$_} / length($seq{$taxid})) * $abundance{$taxid}	} keys %ntrate;
}

##mutates the sequence based on frequency
sub mutate
{
my ($readnt,$qual,$readLength) = @_;
    	my @readnt = split(//, $readnt);
    	my @qual = split(/\t/, $qual);		#stores the probability
    	my $loc = 1;				#this is loc of buffered seqeunce in case of deletion event
    	
    	my @outputsequence;	#store sequence
	
	for(my $i=0; $i < $readLength; $i++) { 
	my $qualityscore = $qual[$i];	
	if(rand()<$qualityscore){
	#MUTATION
			#INDEL EVENT #####################################
			if(rand()<$indelrate)	{
			##################################################	
				if(int(rand(2))){ 	#INSERTION
					my $extra = rand_nt();
					my $extra.= $readnt[$loc];
					push @outputsequence, $extra;
					$loc++
				}else{			#DELETION	
					$loc++;
					push @outputsequence, $readnt[$loc];
					$loc++
				}
			##################################################	
			
			#SUBSTITUTION EVENT ##############################
			}else	{ 						
			##################################################	
			    my $substitutedNT = rand_nt(); 
			    push @outputsequence, $substitutedNT;
					$loc++
			    	} 	
			##################################################	
	}else{
	#NO MUTATION
	    push @outputsequence, $readnt[$loc];
					$loc++;
	    	}
	}
	my $finalseq = join '', @outputsequence;
	return $finalseq;
}

##chooses the random ATGC to change
sub rand_nt
{
	my $rand = rand;
	return 'a' if($rand<$globalnt{'a'});
	return 't' if($rand<$globalnt{'a'}+$globalnt{'t'});
	return 'g' if($rand<$globalnt{'a'}+$globalnt{'t'}+$globalnt{'g'});
	return 'c' if($rand<$globalnt{'a'}+$globalnt{'t'}+$globalnt{'g'}+$globalnt{'c'});
	return 'n';
}

#chooses a random taxa
sub choosetaxa 
{
    my $rand = rand();
    foreach (@taxid) { 
	    if($rand > $range{$_}[0] & $rand <= $range{$_}[1])
	    {
	    return $_;
	    }
	}
}

sub processQual
{
my ($qual) = @_;
	my @qual = map { ord($_) - 33 } split('',$qual);	#convert ASCII to indexNum
	my $qualitystring = join "\t", map {$score{$_}} @qual;
return $qualitystring;
}

sub writeSequence
{
my (	$id, $sequence, $qual, $nameOfSequence, $start, $readLength, 
	$outputfile, $type, $filehandle, $filehandle2, $pair) = @_;
my $header = "simuREAD_$id|gi|$gitaxid{$nameOfSequence}|taxID|$nameOfSequence|loc|$start-";
$header .= $start+$readLength."|output|$outputfile/$pair";
if($type eq 'fastq') 
	{ 
	    if($pair == 1) { 
		say $filehandle '@'.$header;
	    }else{
		say $filehandle2 '@'.$header;
	    }
		say $filehandle $sequence;	#sequence
		say $filehandle "+";
		say $filehandle $qual;
	}else{
	  say $filehandle '>'.$header;
	  say $filehandle $sequence;	#sequence
	}
}

sub processPairReads
{
my ($taxaofchoice,$readLength,$qualONE,$qualTWO,$asciiQualONE, $asciiQualTWO) = @_;
#choose location
	my $insertSize = int(random_normal(150, 5));

	#xc suggested to set min size in case probability kills u 

	my $genomelocation = int(rand($leng{$taxaofchoice}-$insertSize+1));

#read one 
	my $readnt = substr($seq{$taxaofchoice}, $genomelocation,$insertSize);	#problem
	my $newsequence = mutate($readnt, $qualONE,$readLength);   	
#	say "my mutated sequence ", $newsequence;


	my $pair = 1;
	
	my $type ='fasta';
	writeSequence($id,$newsequence,$asciiQualONE,$taxaofchoice,$genomelocation,
	$readLength,$outputfile,$type, $fastaOutput, $fastaOutput, $pair);

	my $type ='fastq';
	writeSequence($id,$newsequence,$asciiQualONE,$taxaofchoice,$genomelocation,
	$readLength,$outputfile,$type, $fastqOutputONE, $fastqOutputTWO, $pair);
	$id++;
#read two
	$readnt = scalar reverse $readnt;
	$newsequence = mutate($readnt, $qualTWO, $readLength);   	

	my $pair = 2;
	
	my $type ='fasta';
	writeSequence($id, $newsequence,$asciiQualTWO,$taxaofchoice,$genomelocation,
	$readLength,$outputfile,$type, $fastaOutput,$fastaOutput, $pair);

	my $type ='fastq';
	writeSequence($id, $newsequence,$asciiQualTWO,$taxaofchoice,$genomelocation,
	$readLength,$outputfile,$type, $fastqOutputONE, $fastqOutputTWO, $pair);
	$id++;
}

sub read_fasta_file {
my ($fasta_file) = @_;
   my (%seq,$taxid, $gi);
	
	open my $in, '<', $fasta_file || die "can't open $fasta_file: $!\n";
while (<$in>) {
   	next if(m/^\s*$/); #remove blank lines
	# read fasta
	if(m/>\s*?(\S+)/){
    		($taxid, $gi) = (split(/\|/, $1))[1,3];	#the taxid of the sequence/genome
		$gitaxid{$taxid} = $gi;
	}else{
		chomp;
		$seq{$taxid} .= $_;
	}
}
close $in;
return \%seq;
}




sub gaussian_rand {
    my ($u1, $u2);  # uniformly distributed random numbers
    my $w;          # variance, then a weight
    my ($g1, $g2);  # gaussian-distributed numbers

    do {
        $u1 = 2 * rand() - 1;
        $u2 = 2 * rand() - 1;
        $w = $u1*$u1 + $u2*$u2;
    } while ( $w >= 1 );

    $w = sqrt( (-2 * log($w))  / $w );
    $g2 = $u1 * $w;
    $g1 = $u2 * $w;
    # return both if wanted, else just one
#    return wantarray ? ($g1, $g2) : $g1;
return $g1;
}

sub random_normal
{
my ($mean, $sd) = @_;
my $num = gaussian_rand();
my $newnum = ($num * $sd) + $mean;
}
