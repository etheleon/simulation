#!/usr/bin/env perl 

use strict;
use warnings;
use lib "/export2/home/uesu/perl5/lib/perl5"; 
use autodie;
use v5.10;
use Bio::SeqIO;


use lib "/export2/home/uesu/perl5/lib/perl5";
die "usage: $0 <genera.names.file>\n" unless $#ARGV==0;
#outsim.0101.missing.txt

use REST::Neo4p;
use REST::Neo4p::Query;
my $cypher = "http://192.168.100.1:7474";
REST::Neo4p->connect($cypher);

my $dbDIR = '/export2/home/uesu/db/refseq/arch2/'.'*';
my $dump = '/export2/home/uesu/db/taxonomy_latest/gi_taxid_nucl.dmp.gz';
my $stmt = <<EOF;
START genus = node:ncbitaxid(taxid={nameoftaxa})
MATCH p=(genus)<-[:childof*0..]-(lowerr:species)<-[:childof*0..]-(lowest)
return 
    genus.taxid as originID, 
    lowest.taxid as targetID
EOF

my @refseq = glob "$dbDIR";
#Hash
my (%missingGenera, %missingChildren, %gihash); 
# missingGenera 	GeneraOfInterest::NULL
# missingChildren	ChildOfInterest::parentGenus

##################################################
#Part 1
say "Storing genera with no complete genoems into hash";
##################################################

open my $input, '<', $ARGV[0];
while(<$input>) { 
    if($. !=1){
    chomp;
    $missingGenera{(split(/\t/))[0]}++;
    }}

##################################################
say "#Genera of interest stored in hash ...\n#Now querying graphDB for children";
##################################################

foreach my $taxaname (keys %missingGenera) {
	my $query = REST::Neo4p::Query->new($stmt,{ nameoftaxa => $taxaname});
	$query->execute;
	while (my $row = $query->fetch){
	    my $genera 	= $row ->[0];
	    my $child 	= $row ->[1];
	    $missingChildren{$child} = $genera;
	}
}

#in case the gi is tagged to genus
foreach (keys %missingGenera) { 
	$missingChildren{$_} = $_;
}
say "Found all ",scalar keys %missingChildren," of their children";

open GITAX, "zcat $dump | ";
while(<GITAX>){
    chomp;
    my($gi, $taxid) = split(/\t/);
    if (exists $missingChildren{$taxid})
    	{	#links GI with parentGenus
    $gihash{$gi} = join("\t" , $taxid, $missingChildren{$taxid});
    	}
}
close GITAX;

say "#Finished stored gis of genera with missing Genomes.\nTotal",scalar keys %gihash;

#loop through all of the refseq files and record the refseq ID
open my $output, ">", 	"out/sim.0102.out";
open my $output2, ">", 	"out/sim.0102.out2";
my  $seqout= Bio::SeqIO->new( -format => 'Fasta', -file => '>out/sim.0102.output.fna');

say $output (join "\t", qw/gi taxid parentGenus refid/);	#HEADER
map {$missingGenera{$_} = 0} keys %missingGenera;		#resets count to zero

foreach my $fna (@refseq){
my $in  = Bio::SeqIO->new(-file => "zcat $fna |", -format => 'Fasta');
		    while (my $seq = $in->next_seq ){
			my $sequenceID=$seq->display_id;
			my ($gi, $refid) = (split(/\|/, $sequenceID))[1,3];
			if (exists $gihash{$gi}) 
			{
			    my $parentGenus = (split(/\t/,$gihash{$gi}))[1];
			    $missingGenera{$parentGenus} += $seq->length;
			    say $output join("\t", $gi,$gihash{$gi},$refid);
			    $seqout->write_seq($seq);
			}
				    }
				}
map {say $output2 "$_\t$missingGenera{$_}" } keys %missingGenera;
