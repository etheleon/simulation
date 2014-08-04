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

my $stmt = <<EOF;
match (genus:genus) 
where genus.name = {nameoftaxa}
with  genus 
match p=(genus)<-[:childof*0..]-(lowerr:species)<-[:childof*0..]-(lowest)
return 
    genus.taxid as originID, 
    genus.name as originName, 
    lowest.taxid as targetID, 
    lowest.name as targetName,
    head(labels(lowest)) as rank
EOF

say (join "\t",qw/originID originName targetID targetName rank/);
foreach my $taxaname (keys %topgenera) {
	my $query = REST::Neo4p::Query->new($stmt,{ nameoftaxa => $taxaname});
	$query->execute;
	while (my $row = $query->fetch) {
	    say (join ("\t", $row->[0], $row->[1], $row->[2], $row->[3], $row->[4]));
	}
}
