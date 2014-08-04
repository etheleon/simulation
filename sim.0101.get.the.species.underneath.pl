#!/usr/bin/env perl

use strict;
use v5.10;
use autodie;
use lib "/export2/home/uesu/perl5/lib/perl5";
die "usage: $0 <genera.names.file>\n" unless $#ARGV==0;
#Aim of this is to parse

use REST::Neo4p;
use REST::Neo4p::Query;
REST::Neo4p->connect("http://192.168.100.1:7474");

my %topgenera;

while(<>) { 
    if($. !=1){ 
    chomp;
    $topgenera{(split(/\t/))[0]}++;
    }}
say "#Genera are all stored in hash ...\n#Now querying graphDB";

#my $stmt='start basetaxa=node:ncbitaxid(taxid={taxids}) match basetaxa-[:childof*]->(genus:`genus`) return genus.taxid';

#First match
my $stmt='match (n:genus) where n.name ={nameoftaxa} return n.taxid';
foreach my $taxaname (keys %topgenera) {
	my $query = REST::Neo4p::Query->new($stmt,{ nameoftaxa => $taxaname});
	$query->execute;
	while (my $row = $query->fetch) {
		$topgenera{$taxaname} = $row->[0];
	}
}
say '#Assigned taxid to genera';

#2nd match find all taxa at and below species associated with the genus of interest
my $stmt2 = <<EOF;
START basetaxa=node:ncbitaxid(taxid={taxids}) 
MATCH basetaxa<-[:childof*0..]-(lower)<-[:childof*0..]-(lowest)
WHERE lower:species 
RETURN labels(lowest), lowest.taxid, lowest.name
EOF

say "taxid\trank\t\tparentgenera\tchildname";	#header

foreach my $parenttaxid (keys %topgenera){
	my $query = REST::Neo4p::Query->new($stmt2,{ taxids => $topgenera{$parenttaxid}});
	$query->execute;
	while (my $row = $query->fetch) {
		my ($rank, $taxid, $name);
#	   #extract the rank##################################################
	    	foreach (@{$row->[0]}) {		#this column has multiple items hence its an array
	    	    if (!/Taxon/){			#there's labels called "Taxon"
	    	    	$rank = $_;
	    	    }}
#	   ###################################################################
		$taxid = $row->[1];
		$name  = $row->[2];
		say "$taxid\t$rank\t$topgenera{$parenttaxid}\t$name";
	}
}
