#!/usr/bin/env perl 

use strict;
use warnings;
use autodie;
use v5.20;

use lib "/export2/home/uesu/perl5/lib/perl5";
die "usage: $0 <genera.names.file>\n" unless $#ARGV==0;

use REST::Neo4p;
use REST::Neo4p::Query;

##################################################
#Init
##################################################
#my $dbDIR = '/export2/home/uesu/db/refseq/arch2/'.'*';
my $dbDIR = '/Volumes/altrepo/water/db/'.'*';
#my $dump = '/export2/home/uesu/db/taxonomy_latest/gi_taxid_nucl.dmp.gz'
my $dump = '/Users/uesu/gi_taxid_nucl.dmp.gz';
my %gihash;
my $stmt = <<EOF;
START genus = node:ncbitaxid(taxid={nameoftaxa})
MATCH p=(genus)<-[:childof*0..]-(lowerr:species)<-[:childof*0..]-(lowest)
return 
    genus.taxid as originID, 
    lowest.taxid as targetID, 
EOF
my $cypher = "http://localhost:7474";
my @refseq = <$dbDIR>;

#Hash
my (%missingGenera, %missingChildren); 
# missingGenera 	GeneraOfInterest::
# missingChildren	ChildOfInterest::gi

#Misc
REST::Neo4p->connect($cypher);

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

open GITAX, "zcat $dump | ";
while(<GITAX>){
    chomp;
    my($gi, $taxid) = split(/\t/);
    $gihash{$gi} = $taxid if (exists $missingGenera{$taxid});	#links GI with parentGenus
}
close GITAX;
say "#Finished stored gis of genera with missing Genomes.\nTotal",scalar keys %gihash;

#loop through all of the refseq files and record the refseq ID
say "gi\ttaxid\trefid";	#HEADER
open my $output, ">", "out/sim.0105.out";
foreach my $fna (@refseq){
		    open my $fasta, '<', $fna; 
		    while(<$fasta>)
		    		{ 
			my ($gi, $refid) = (split(/\|/)[1,3] if /^\>/;
			if(exists $gihash{$gi}) { say "$gi\t$gihash{$gi}\t$refid" } 
				}
		}
#my $numberofjobs = scalar @refseq;

#taxID	giID	refseqID
#
