#!/export2/home/uesu/local/perl52/bin/perl

use v5.20;
use feature 'signatures';
no warnings 'experimental';
use autodie;
use Bio::SeqIO;
use lib "/export2/home/uesu/perl5/lib/perl5";

die "usage: $0 chosenGenomes.fna template.fq output indel-rate abundanceInfo
		(indel-rate: proportion of indel over all errors; default=0.1)
		eg. $0 out/sim.0200.out.fna /export2/ulu_pandan2011/data/batch1_gDNA_illumina/filtered_fastq/s_1_1.filtered.fastq out/sim.0300.output 0.01 \n" unless $#ARGV>=2;

##################################################
# Part0: Init
say "Initializing ...";
##################################################

#Hashes###
my (%seq,%leng,%score,%ntrate,%abundance,%range);
#	%seq 	taxaID::concatenated sequence of chromosomes
#	%leng 	taxaID::length of the sequence
#	%ntrate	nucleotide::the frequencies
##########

#Arrays###
my @taxid;						#stores taxid in decreasing abundances
##########

#Scalar##
my ($totalleng,$totalAbu);				#length of 
my $id = 0;						#For adding to fastq header
my $indelrate=0.01;				 	#the INDEl error rate:percentage of errors to be indel
my $outputfile = $ARGV[2]; $outputfile =~ s/.+\///; 	#the outputfile name
##########

#Phred error probabilities
map {$score{$_}=phred($_)} 0..100; 			#Probability 


##################################################
# Part1: read genomes 
say "Reading Genomes ...";
##################################################

my $genome = Bio::SeqIO->new(-file => "$ARGV[0]", -format => 'Fasta');
while (my $sequence = $genome->next_seq ){
	my $ntsequence 	=$sequence->seq;			#nucleotide sequence 
    	my $taxid = (split(/\|/, $sequence->display_id))[1];	#the taxid of the sequence/genome
	$seq{$taxid}	.=$ntsequence; 				#concatenates sequences of the same taxa tog
	$leng{$taxid}	+=$sequence->length; 			#the length of the concatenated sequences
}$genome->close;

my @chrom = sort keys %seq;

#	$ntrate{'a'} += ($ntsequence=~tr/aA/AA/);	#count the number of ATCGs
#	$ntrate{'t'} += ($ntsequence=~tr/tT/TT/);
#	$ntrate{'g'} += ($ntsequence=~tr/gG/GG/);
#	$ntrate{'c'} += ($ntsequence=~tr/cC/CC/);
#$ntrate{$_}/=$totalleng foreach(keys %ntrate);		#calculate nucleotide frequency
#say scalar @chrom." chromosomes with $totalleng base pairs were read."; 	
#map {say "\t$_: $ntrate{$_}"} keys %ntrate;

###################################################
# Part2: abundance
say "reading abundances";
##################################################

open my $abundance, '<', "../out/sim.0201.out.txt";
while(<$abundance>){unless ($. == 1){ 
	my ($abu, $taxid) = (split(/\t/))[1,3];
	$abundance{$taxid} = $abu
}}

#abundance summary
map { $totalAbu += $abundance{$_} } keys %abundance;	#total abundance 
map { $abundance{$_} = ($abundance{$_}/$totalAbu)} keys %abundance;

#taxid - range creation:: TAXID: (lower, upper)
@taxid= sort { $abundance{$b} <=> $abundance{$a} } keys %abundance;
my @zipped = map {($taxid[$_], $abundance{$taxid[$_]})} 0 .. $#taxid;

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

#foreach (@taxid) { say "$_\t$range{$_}[0]\t$range{$_}[1]" }

##################################################
# reading quality scores
say "Reading quality scores & simulating shotgun ...";
##################################################

open my $output, '>', $ARGV[2];				#OUTPUT
open my $fastq, '<', $ARGV[1];
while(my $header = <$fastq>){
    	chomp $header;
	<$fastq>;<$fastq>;
	
	$id++;							#the nth fastq sequence processed

#quality
	my $qual =  <$fastq>;
	chomp $qual;
	my @qual = map { ord($_) - 33 } split('',$qual);	#stores quality scores into array
	my $leng = @qual;					#the length of the fastq read 

#Choosing taxa and position
	my $taxaofchoice = choosetaxa(@taxid);
	say $taxaofchoice;
	my $readsize 		=  	$leng * 2;			#nt size be plucked from genome
	#choose loc
	my $genomelocation 	= 	int(rand($seq{$taxaofchoice}-$readsize+1));
	$genomelocation = int(rand($taxaofchoice-$readsize+1)) while(! fitChromosome($genomelocation, $readsize));	#reassign 

	my $source = fitChromosome($genomelocation, $readsize);
	my @source = @$source;
	# read the sequence from the fitted genomic region
	my @readnt = split('',substr($seq{$source[0]}, $source[1],$readsize));
	
#	# create a new shotgun sequencing read
	my $newread;
#	my $index=0;
	my $newsequence = map { mutate(@readnt,@qual, $_)} 0..$#qual;   	#generates new sequence
	
	say $output "";			#header
	say $output $newsequence;	#sequence
}
close $fastq;


#open QUAL,'<', $ARGV[1];
#while(my $h = <QUAL>)
#{
#	chomp $h;
#	my $seq = <QUAL>;
#	my $h2 = <QUAL>;
#	my $qual = <QUAL>;
#	chomp $qual;	

#	my @qual = map { ord($_) - 33} split('', $qual);	#converts the ASCII chars into numbers 
#	my $leng = @qual;					#the length of the fastQ read 
#	next unless $leng;
#
#	$id++;
#	# find a starting point for a sequencing read; 
#	# it must fit into the chromosome, or the length between 
#	# the point and end of chromosome must be longer than the read 
#	my $size = $leng + int($leng*$indelrate);			#the size of the read
#	my $genomelocation = int(rand($totalleng-$size+1));		#location of the concatenated genomes
#
#	$genomelocation = int(rand($totalleng-$size+1)) while(! fitChromosome($genomelocation, $size));
#	my $source = fitChromosome($genomelocation, $size);
#	my @source = @$source;
#	# read the sequence from the fitted genomic region
#	my @readnt = split('',substr($seq{$source[0]}, $source[1],$size));
#	
#	# create a new shotgun sequencing read
#	my $newread;
#	my $index=0;
#	# simulate nucleotide by nucleotide
#	foreach(my $i=0; $i<=$#qual; $i++)
#	{
#		# error probability for this site
#		my $qual = $score{$qual[$i]};
#		# if random probablity fall in the error probality
#		if(rand()<$qual)
#		{

##If you do not want to include indel estimation #################################################
## if random chance says it should be a indel
#			if(rand()<$indelrate)
#			{
#				# insertion or deletion
#				my $insert = int(rand(2));
#				if($insert)
#				{
#					# insert a random nt
#					$newread .= rand_nt();
#					next;
#				}else
#				{
#					# delete one nt by increasing the index
#					$index++;
#				}
#			# if not indel
#			}else
#			{
#####################################################################################################
#				$readnt[$index] = rand_nt();
#			}	#should comment this line out if you want to remove indel simulation
#		}
#		$newread .= $readnt[$index];
#		$index++;
#	}
#
#	say $output '@'."simu_${id}:$source[0]:$source[1]-",($source[1]+$leng),":$outputfile";
#	say $output "$newread\n$h2$qual";
#}		
#close $output;
#say "all done. please check $ARGV[2] for results";
#
####################################################################################################
#Functions
####################################################################################################

sub fitChromosome($loc, $size)
{
	my @location;
	my $leng;
	foreach(@chrom) 	#@chrom stores the number of sequences in the fasta file; loops through from loc 0 till end
	{
		my $locinthis = $loc-$leng;		#the randomlly chosen global loc - leng
		my $thisleng = $leng{$_};	#the length of this chromosome
		$leng += $thisleng;			#stores the 
		next if($locinthis>$thisleng-$size);	#if it fails ie. falls outside of the xsome's length
		@location = ($_,$locinthis,$locinthis+$size-1);	#the sequence, the start location, the end location
		return \@location;
	}
	return;
}

#Calculates error-probability of mutation
sub phred($score)
{
	return 10**($score/(-10));
}

sub mutate(@readnt,@qual,$loc){
	my $qualityscore = $score{$qual[$loc]};	# error probability for this site
	if(rand()<$qualityscore){
			if(rand()<$indelrate){						#INDEL EVENT
				if(int(rand(2))){ 	#INSERTION
					my $extra = rand_nt();
					my $extra.= $readnt[$loc];
					return $extra;
				}else{			#DELETION	
					return "";
				}
			}else{ my $substitutedNT = rand_nt(); return $substitutedNT} 	#SUBSTITUTION EVENT
	}else{return $readnt{$loc}}
}


#	#chooses the random ATGC to change
sub rand_nt
{
	my $rand = rand;
	return 'a' if($rand<$ntrate{'a'});
	return 't' if($rand<$ntrate{'a'}+$ntrate{'t'});
	return 'g' if($rand<$ntrate{'a'}+$ntrate{'t'}+$ntrate{'g'});
	return 'c' if($rand<$ntrate{'a'}+$ntrate{'t'}+$ntrate{'g'}+$ntrate{'c'});
	return 'n';
}

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
