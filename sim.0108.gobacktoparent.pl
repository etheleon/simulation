#!/usr/bin/env perl

use strict;
use v5.10;
use autodie;
use lib "/export2/home/uesu/perl5/lib/perl5";
#Aim of this is to parse

use REST::Neo4p;
use REST::Neo4p::Query;
REST::Neo4p->connect("http://192.168.100.1:7474");

my %taxhash;
open IN, "$ARGV[0]";
while(<IN>){ 
    chomp;
    unless ($. == 1){
    my $taxid = (split(/\s+/))[1];
    $taxhash{$taxid}++;
    }
}
close IN;

my $stmt = <<EOF;
START basetaxa=node:ncbitaxid(taxid={taxids}) 
MATCH basetaxa-[:childof*0..]->(lower)
WHERE lower:genus
RETURN distinct lower.taxid
EOF

my %speciesgenus;
foreach my $taxaname (keys %taxhash) {
	my $query = REST::Neo4p::Query->new($stmt,{ taxids => $taxaname});
	$query->execute;
	while (my $row = $query->fetch) {
	    my $genus = $row->[0];
	    $speciesgenus{$taxaname} = $genus;
	}
}
#foreach (keys %speciesgenus) { 
#    say "$_\t$speciesgenus{$_}"
#}

open OUT, ">",$ARGV[1];
open IN, $ARGV[0];
while(<IN>){ 
    chomp;
    unless ($. == 1){
    my $taxid = (split(/\s+/))[1];
    say OUT "$_"."\t"."$speciesgenus{$taxid}";
    }
}
close IN;
