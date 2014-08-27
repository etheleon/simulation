#!/usr/bin/env perl

use strict;
use v5.10;
use autodie;
use lib "/export2/home/uesu/perl5/lib/perl5";
die "usage: $0 <genera.names.file>\n" unless $#ARGV==0;	#eg. sim.0101.out2.txt
#Aim of this is to parse

use REST::Neo4p;
use REST::Neo4p::Query;
REST::Neo4p->connect("http://192.168.100.1:7474"); #need to update the internal DB

my %topgenera;

while(<>) { 
    if($. !=1){ 
    chomp;
    $topgenera{(split(/\t/))[1]}++;
    }}

open my $output, ">", "out/sim.0200.out.txt";
say $output "#Genera are all stored in hash ...\n#Now querying graphDB";
#foreach (keys %topgenera) { say };

my $stmt = <<EOF;
START genus = node:ncbitaxid(taxid={nameoftaxa})
MATCH p=(genus)<-[:childof*0..]-(lowerr:species)<-[:childof*0..]-(lowest)
return 
    genus.taxid as originID, 
    genus.name as originName, 
    lowest.taxid as targetID, 
    lowest.name as targetName,
    head(labels(lowest)) as rank
EOF

say $output (join "\t",qw/originID originName targetID targetName rank/);
foreach my $taxaname (keys %topgenera) {
	my $query = REST::Neo4p::Query->new($stmt,{ nameoftaxa => $taxaname});
	$query->execute;
	while (my $row = $query->fetch){
	    say $output (join ("\t", $row->[0], $row->[1], $row->[2], $row->[3], $row->[4]));
	}
}
