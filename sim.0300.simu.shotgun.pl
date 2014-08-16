#!/export2/home/uesu/local/perl52/bin/perl

use v5.20;
use feature 'signatures';no warnings 'experimental';
use autodie;
use Bio::SeqIO;
use lib "/export2/home/uesu/perl5/lib/perl5";

die "usage: $0 chosenGenomes.fna template.fq output indel-rate abundanceInfo
		(indel-rate: proportion of indel over all errors; default=0.1)\n" unless $#ARGV>=2;

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
my (@taxid,@zipped);						#stores taxid in decreasing abundances
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
say "2. Reading Genomes ...";
##################################################

my $genome = Bio::SeqIO->new(-file => "$ARGV[0]", -format => 'Fasta');
while (my $sequence = $genome->next_seq ){
	my $ntsequence 	=$sequence->seq;			#nucleotide sequence 
    	my $taxid = (split(/\|/, $sequence->display_id))[1];	#the taxid of the sequence/genome
	$seq{$taxid}	.=$ntsequence; 				#concatenates sequences of the same taxa tog
}

###################################################
# Part2: abundance
say "3. Reading abundances.. ";
##################################################

open my $abundance, '<', "../out/sim.0201.out.txt";
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
# reading quality scores
say "4. Reading quality scores & simulating shotgun ...";
##################################################

open my $output, 	'>', 	$ARGV[2];				#OUTPUT
open my $fastq, 	'<', 	$ARGV[1];
while(my $header = <$fastq>){
    	chomp $header;
	<$fastq>;<$fastq>;
	$id++;							#the nth fastq sequence processed

#quality
	my $qual =  <$fastq>;
	chomp $qual;
	my @qual = map { ord($_) - 33 } split('',$qual);	#stores quality scores into array
	my $qualitystring = join "\t", @qual;
	my $leng = @qual;					#the length of the fastq read 

#Choosing taxa and loc to pluck sequence out from
	my $taxaofchoice = choosetaxa(@taxid);
	say "\t Chose taxaID: $taxaofchoice";
	my $readsize 		=  	$leng * 2;			#nt size be plucked from genome

	#choose loc
	say "Length of sequence:",length $seq{$taxaofchoice};
	my $genomelocation 	= 	int(rand(length($seq{$taxaofchoice})-$readsize+1));

	$genomelocation = int(rand($taxaofchoice-$readsize+1)) while(! fitChromosome($taxaofchoice,$genomelocation, $readsize));	#reassign if it doesnt fit
	say "Location to sample: $genomelocation";

	my $source = fitChromosome($taxaofchoice, $genomelocation, $readsize);
	my @source = @$source;
	
	# read the sequence from the fitted genomic region
	my $readnt = substr($seq{$source[0]}, $source[1],$readsize);	#problem
	    say "\tthe string's length is: ",length $readnt;

	# create a new shotgun sequencing read
	my $newread;
	my $index=0;
	my $newsequence = join '', (map { mutate($readnt,$qualitystring,$_)} 0..$#qual);   	#generates new sequence
	
	say $output '>'."simu_${id}:$source[0]:$source[1]-",($source[1]+$leng),":$outputfile";
#	eg >simu_taxid|111|GI|0010101|ref|NC_00001:0-10:outputfilename
	say $output $newsequence;	#sequence
}
close $fastq;


####################################################################################################
#Functions
####################################################################################################

#returns a suitable location
sub fitChromosome($taxid, $loc, $size)	
{
	my @location;
	if( (length($seq{$taxid}) - $loc) > $size)
	{
	    @location = ($taxid, $loc, $loc+$size-1);
	    return \@location;
	}
	return;
}



#Calculates error-probability of mutation
sub phred($score)
{
	return 10**($score/(-10));
}

sub ntfreq($taxid)
{
    	my %ntrate;
	$ntrate{'a'} += ($seq{$taxid}=~tr/aA/AA/);	#count the number of ATCGs
	$ntrate{'t'} += ($seq{$taxid}=~tr/tT/TT/);
	$ntrate{'g'} += ($seq{$taxid}=~tr/gG/GG/);
	$ntrate{'c'} += ($seq{$taxid}=~tr/cC/CC/);
	map {	$globalnt{$_} += ($ntrate{$_} / length $seq{$taxid}) * $abundance{$taxid}	} keys %ntrate;
}

#mutates the sequence based on frequency
sub mutate($readnt,$qual,$loc){
    	my @readnt = split(//, $readnt);
    	my @qual = split(/\t/, $qual);
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
	}else{return $readnt[$loc]}
}


#	#chooses the random ATGC to change
sub rand_nt
{
	my $rand = rand;
	return 'a' if($rand<$globalnt{'a'});
	return 't' if($rand<$globalnt{'a'}+$globalnt{'t'});
	return 'g' if($rand<$globalnt{'a'}+$globalnt{'t'}+$globalnt{'g'});
	return 'c' if($rand<$globalnt{'a'}+$globalnt{'t'}+$globalnt{'g'}+$globalnt{'c'});
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

say "All done. Please check $ARGV[2] for results";
