#!/usr/bin/perl -w 

use 5.012;
use autodie;
use lib "/export2/home/uesu/perl5/lib/perl5";

die "usage: $0 template.fasta shotgun.phred.qual out.file indel-rate abundanceinformation
\t\t(indel-rate: proportion of indel over all errors; default=0.1)\n" unless $#ARGV>=2;

my %leng;
my %seq;
my %score;
my %ntrate;

my $chrom;
my $indelrate=0.01;
my $totalleng;
my $id = 0;
my $file = $ARGV[2];
$file =~ s/.+\///;

#OUTPUT
open OUT, '>', $ARGV[2];

# indel rate: percentage of errors to be indel
$indelrate = $ARGV[3] if $#ARGV>=3;

# pre-calculate error-probability based on phred score
$score{$_}=phred($_) foreach(0..100);

# read the genome fasta file
say "reading genome file...";
#fasta
open my $genome, '<', $ARGV[0];
while(<$genome>)
{
	next if(m/^\s*$/); 	# skip emputy lines
	if(m/>\s*?(\S+)/)	# read fasta 
	{
		$chrom = $1;	#header
	}else
	{
		$seq{$chrom} .= $_;	#concatenates the sequences following the header into a single string
	}
}

# calculate nucleotide frequency and length of the genome
say "analyzing genome..."; 

foreach my $chrom (keys %seq)	#foreach sequence
{
	$seq{$chrom} =~ s/\s//g;			#remove space chars
	$leng{$chrom} = length $seq{$chrom};		#store the length of the sequence
	$totalleng += $leng{$chrom};			#stores the total length of all the genomes into $totalength
	$ntrate{'a'} += ($seq{$chrom}=~tr/aA/AA/);	#count the number of ATCGs
	$ntrate{'t'} += ($seq{$chrom}=~tr/tT/TT/);
	$ntrate{'g'} += ($seq{$chrom}=~tr/gG/GG/);
	$ntrate{'c'} += ($seq{$chrom}=~tr/cC/CC/);
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
	my $leng = @qual;
	next unless $leng;

	$id++;							#the nth FASTQ sequence

	# find a starting point for a sequencing read; 
	# it must fit into the chromosome, or the length between 
	# the point and end of chromosome must be longer than the read 
	my $size = $leng + int($leng*$indelrate);
	my $genomelocation = int(rand($totalleng-$size+1));
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
	say OUT '@'."simu_${id}:$source[0]:$source[1]-",($source[1]+$leng),":$file";
	say OUT "$newread\n$h2$qual";
}		
close OUT;
say "all done. please check $ARGV[2] for results";

#Functionas
sub fitChromosome
{
	my $loc = shift;
	my $size = shift;
	my @location;
	my $leng;
	foreach(0..$#chrom)
	{
		my $locinthis = $loc-$leng;
		my $thisleng = $leng{$chrom[$_]};
		$leng += $thisleng;
		next if($locinthis>$thisleng-$size);
		@location = ($chrom[$_],$locinthis,$locinthis+$size-1);
		return \@location;
	}
	return;
}

sub phred
{
	#QV = - 10 * log_10( P_e )
	my $score = shift;
	return 10**($score/(-10));
}

sub rand_nt
{
	my $rand = rand;
	return 'a' if($rand<$ntrate{'a'});
	return 't' if($rand<$ntrate{'a'}+$ntrate{'t'});
	return 'g' if($rand<$ntrate{'a'}+$ntrate{'t'}+$ntrate{'g'});
	return 'c' if($rand<$ntrate{'a'}+$ntrate{'t'}+$ntrate{'g'}+$ntrate{'c'});
	return 'n';
}
