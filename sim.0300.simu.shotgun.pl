#!/export2/home/uesu/local/perl52/bin/perl
#testing
use v5.20;
use feature 'signatures';
no warnings 'experimental';
use autodie;
use Bio::SeqIO;
use lib "/export2/home/uesu/perl5/lib/perl5";

die "usage: $0 chosenGenomes.fna template.fq output indel-rate abundanceInfo
\t\t(indel-rate: proportion of indel over all errors; default=0.1)
\t\teg. $0 out/sim.0200.out.fna\n /export2/ulu_pandan2011/data/batch1_gDNA_illumina/filtered_fastq/s_1_1.filtered.fastq out/sim.0300.output 0.01 \n" unless $#ARGV>=2;

##################################################
#Init
##################################################
my %leng;
my %seq;
my %score;
# pre-calculate error-probability based on phred score
$score{$_}=phred($_) foreach(0..100);

my %ntrate;

my $chrom;
my $indelrate=0.01;
my $totalleng;
my $id = 0;
my $file = $ARGV[2];
$file =~ s/.+\///;

#OUTPUT
open my $output, '>', $ARGV[2];

# indel rate: percentage of errors to be indel
$indelrate = $ARGV[3] if $#ARGV>=3;

##################################################
# Part1: read genomes 
say "Reading Genomes ...";
##################################################

my $genome = Bio::SeqIO->new(-file => "$ARGV[0]", -format => 'Fasta');
while (my $seq = $genome->next_seq ){


}

#open my $genome, '<', $ARGV[0];
#while(<$genome>)
#{
	next if(m/^\s*$/); 	# skip emputy lines
	if(m/>\s*?(\S+)/)	# read fasta 
	{$chrom = $1}else{$seq{$chrom} .= $_}
}


##################################################
# Part2: read abundances
say "analyzing genome..."; 
##################################################

foreach my $genome (keys %seq)	#foreach sequence
{
	$seq{$genome} =~ s/\s//g;			#remove space chars
	$leng{$genome} = length $seq{$genome};		#store the length of the sequence
	$totalleng += $leng{$genome};			#stores the total length of all the genomes into $totalength
	$ntrate{'a'} += ($seq{$genome}=~tr/aA/AA/);	#count the number of ATCGs
	$ntrate{'t'} += ($seq{$genome}=~tr/tT/TT/);
	$ntrate{'g'} += ($seq{$genome}=~tr/gG/GG/);
	$ntrate{'c'} += ($seq{$genome}=~tr/cC/CC/);
}

$ntrate{$_}/=$totalleng foreach(keys %ntrate);		#calculate nucleotide frequency
my @chrom = sort keys %seq;				#stores the keys into array	#### I need to mod this to take into accound the relative abundance

say @chrom." chromosomes with $totalleng base pairs were read."; 	
say "\t$_: $ntrate{$_}" foreach(keys %ntrate);

# reading quality scores
say "reading quality scores & simulating shotgun ...";

#sequencing quality files
open QUAL,'<', $ARGV[1];
while(my $h = <QUAL>)
{
	chomp $h;
	my $seq = <QUAL>;
	my $h2 = <QUAL>;
	my $qual = <QUAL>;
	chomp $qual;	

	my @qual = map { ord($_) - 33} split('', $qual);	#converts the ASCII chars into numbers 
	my $leng = @qual;					#the length of the fastQ read 
	next unless $leng;

	$id++;							#the nth FASTQ sequence

	# find a starting point for a sequencing read; 
	# it must fit into the chromosome, or the length between 
	# the point and end of chromosome must be longer than the read 
	my $size = $leng + int($leng*$indelrate);			#the size of the read
	my $genomelocation = int(rand($totalleng-$size+1));		#location of the concatenated genomes
	$genomelocation = int(rand($totalleng-$size+1)) while(! fitChromosome($genomelocation, $size));
	my $source = fitChromosome($genomelocation, $size);
	my @source = @$source;
	# read the sequence from the fitted genomic region
	my @readnt = split('',substr($seq{$source[0]}, $source[1],$size));
	
	# create a new shotgun sequencing read
	my $newread;
	my $index=0;
	# simulate nucleotide by nucleotide
	foreach(my $i=0; $i<=$#qual; $i++)
	{
		# error probability for this site
		my $qual = $score{$qual[$i]};
		# if random probablity fall in the error probality
		if(rand()<$qual)
		{

#If you do not want to include indel estimation #################################################
# if random chance says it should be a indel
			if(rand()<$indelrate)
			{
				# insertion or deletion
				my $insert = int(rand(2));
				if($insert)
				{
					# insert a random nt
					$newread .= rand_nt();
					next;
				}else
				{
					# delete one nt by increasing the index
					$index++;
				}
			# if not indel
			}else
			{
####################################################################################################
				$readnt[$index] = rand_nt();
			}	#should comment this line out if you want to remove indel simulation
		}
		$newread .= $readnt[$index];
		$index++;
	}
	say $output '@'."simu_${id}:$source[0]:$source[1]-",($source[1]+$leng),":$file";
	say $output "$newread\n$h2$qual";
}		
close $output;
say "all done. please check $ARGV[2] for results";

#Functions
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
		@location = ($_,$locinthis,$locinthis+$size-1);	#the length of the chromosome, the start location, the end location
		return \@location;
	}
	return;
}

#Calculates probability of mutation
sub phred($score)
{
	return 10**($score/(-10));
}

#chooses the random ATGC to change
sub rand_nt
{
	my $rand = rand;
	return 'a' if($rand<$ntrate{'a'});
	return 't' if($rand<$ntrate{'a'}+$ntrate{'t'});
	return 'g' if($rand<$ntrate{'a'}+$ntrate{'t'}+$ntrate{'g'});
	return 'c' if($rand<$ntrate{'a'}+$ntrate{'t'}+$ntrate{'g'}+$ntrate{'c'});
	return 'n';
}
