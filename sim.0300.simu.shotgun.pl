#!/usr/bin/env perl

use v5.10;
use lib "/export2/home/uesu/perl5/lib/perl5";
use autodie;
use Bio::SeqIO;
use Math::Random;

die "usage: $0 chosenGenomes.fna template.fq output abundanceInfo 
		(indel-rate: proportion of indel over all errors; default=0.1)\n" unless $#ARGV>=3;

##################################################
# Part0: Init
say "1. Initializing ...";
##################################################

#Hashes###
my (%seq,%leng,%score,%ntrate,%abundance,%range,%globalnt);
#	%seq 	taxaID::concatenated sequence of chromosomes
#	%leng 	taxaID::length of the sequence
#	%ntrate	nucleotide::the frequencies
#	%abundance taxaID::abundance
#	%range 	taxaID::the lower and upper bounds
#	%globalnt nucleotide::frequency
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
my $fastqfile  = $ARGV[1]; $fastqfile =~ s/.+\///; 	#the fastqfile 
##########

#Phred error probabilities
map {$score{$_}=phred($_)} 0..100; 			#Probability 
##################################################
# Part1: read genomes 
say "2. Reading Genomes ...";
##################################################

my $genome = Bio::SeqIO->new(-file => "$ARGV[0]", -format => 'Fasta');
while (my $sequence = $genome->next_seq ){
	my $ntsequence 	=$sequence->seq;			#nucleotide sequence 
    	my $taxid = join '_', (split(/\|/, $sequence->display_id))[1,3];	#the taxid of the sequence/genome
	$seq{$taxid}	.=$ntsequence; 				#concatenates sequences of the same taxa tog
}

###################################################
# Part2: abundance
say "3. Reading abundances.. ";
##################################################

open my $abundance, '<', $ARGV[3]; #out/sim.0201.out.txt
while(<$abundance>){unless ($. == 1){ 
	my ($abu, $taxid) = (split(/\t/))[1,3];
	$abundance{$taxid} = $abu
}}

#abundance summary
map { $totalAbu += $abundance{$_} } keys %abundance;	#total reads per million abundance 
map { $abundance{$_} = ($abundance{$_}/$totalAbu)} keys %abundance;

#calculate Nucleotide Frequency
map { ntfreq($_) } keys %abundance;

say "\tNucleotide frequencies";
map {say "\t\t$_: $globalnt{$_}"} keys %globalnt;

##taxid - range creation:: TAXID: (lower, upper)
@taxid= sort { $abundance{$b} <=> $abundance{$a} } keys %abundance;	#sort taxid by their abundance
@zipped = map {($taxid[$_], $abundance{$taxid[$_]})} 0 .. $#taxid;

#THINGS TO DO to convert into a mappable function
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

##################################################
# Simulate Shotgun sequence
say "4. Reading quality scores & simulating ...";
##################################################

open my $fastaOutput, 	'>', 	"$ARGV[2]".'.fna';				#FASTA OUTPUT
open my $fastqOutput, 	'>', 	"$ARGV[2]".'.fq';				#FASTQ OUTPUT

open my $fastq1, 	'<', 	$fastqfile."_1.filtered.fastq";
open my $fastq2, 	'<', 	$fastqfile."_2.filtered.fastq";

while (!eof($fastq1) and !eof($fastq2)) {
	<$fastq1>;<$fastq1>;	#sequence, h2
	<$fastq2>;<$fastq2>;	#sequence, h2
#Process Quality
	
	my $qual1 =  <$fastq1>;chomp $qual1;
	my $qual2 =  <$fastq1>;chomp $qual2;

	my $readlength = length($qual);	#template length
	
	my $qualitystring1 = processQual($qual1);
	my $qualitystring2 = processQual($qual2);
	
#Choosing taxa and loc to pluck sequence out from
	my $taxaofchoice = choosetaxa(@taxid);


	processPairReads($taxaofchoice, $readlength, $qualitystring1, $qualitystring2);
}

say "All done. Please check $ARGV[2] for results";
####################################################################################################
#Functions
####################################################################################################

#returns a suitable location
sub fitChromosome
{
my ($taxid, $loc, $size) = @_;
	if( (length($seq{$taxid}) - $loc) > $size)
	{
	    my @location = ($taxid, $loc, $loc+$size-1);
	    return \@location;
	}
	return;
}

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

#mutates the sequence based on frequency
sub mutate
{
my ($readnt,$qual) = @_;
    	my @readnt = split(//, $readnt);
    	my @qual = split(/\t/, $qual);
    	my $loc = 1;	#this is loc of buffered seqeunce in case of deletion event
    	my $i = scalar @qual;		#this is length of qual 
    	
    	my @outputsequence;	#store sequence
	
	for($i; $i > 0; $i--) { 
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

#chooses the random ATGC to change
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
my ($sequence, $qual, $nameOfSequence, $start, $readLength, $outputfile, $type, $filehandle,$pair) = @_;
	if($type eq 'fastq') 
	{ 
	say $filehandle '@'."simu|taxID_gi|$source[0]|loc|$source[1]-",($source[1]+$readLeng),"|output|$outputfile/$pair";
	say $filehandle $newsequence;	#sequence
	say $filehandle "+\n$qual";
	}else
	{
	say $filehandle '>'."simu|taxID_gi|$source[0]|loc|$source[1]-",($source[1]+$readLeng),"|output|$outputfile/$pair";
	say $filehandle $newsequence;	#sequence
	}
}

sub processPairReads
{
my ($taxaofchoice,$readlength,$qual1,$qual2) = @_;
#choose location
	my $insertSize = random_normal(1, 150, 5);
	my $genomelocation 	= 	int(rand(length($seq{$taxaofchoice})-$insertSize+1));
	$genomelocation = int(rand(length($seq{$taxaofchoice})-$insertSize+1)) while(! fitChromosome($taxaofchoice,$genomelocation, $insertSize));	#reassign if it doesnt fit
	my $source = fitChromosome($taxaofchoice, $genomelocation, $insertSize);
	my @source = @$source;

#read one 
	my $readnt = substr($seq{$source[0]}, $source[1],$insertSize);	#problem
	my $newsequence = mutate($readnt, $qual1);   	

	my $pair = 1;
	
	my $type ='fasta';
	writeSequence($newsequence,$qual,$source[0],$source[1],
	$readLength,$outputfile,$type, $fastaOutput, $pair);

	my $type ='fastq';

	writeSequence($newsequence,$qual,$source[0],$source[1],
	$readLength,$outputfile,$type, $fastaOutput, $pair);

#read two
	$readnt = scalar reverse $readnt;
	$newsequence = mutate($readnt, $qual2);   	

	my $pair = 2;
	
	my $type ='fasta';
	writeSequence($newsequence,$qual,$source[0],$source[1],
	$readLength,$outputfile,$type, $fastaOutput, $pair);

	my $type ='fastq';

	writeSequence($newsequence,$qual,$source[0],$source[1],
	$readLength,$outputfile,$type, $fastaOutput, $pair);
}
