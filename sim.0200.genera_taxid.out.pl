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
say "genera.name\ttaxid";
foreach (keys %topgenera){ 
    say "$_\t$topgenera{$_}";
}
