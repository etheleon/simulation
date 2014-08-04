#!/usr/bin/env perl 

use lib "/export2/home/uesu/perl5/lib/perl5";
use strict;
use v5.10;
use autodie;
use Bio::SeqIO;
use IO::Zlib;

#Stores GI-taxid %gihash 
my %gihash;
open GITAX, "zcat /export2/home/uesu/db/taxonomy_latest/gi_taxid_nucl.dmp.gz | ";
while(<GITAX>){
    chomp;
    my($gi, $taxid) = split(/\t/);
	$gihash{$gi}=$taxid;
}
close GITAX;

#Stores refseq GIs of the top100 taxa into %gihash2 (hash of arrays)

#410600 fasta entries
my %gihash2;

open IN, "zcat out/sim.0106.input.gz | ";
open OUT, ">out/sim.0107.out.txt";

say OUT "gi\ttaxid";	#header
while(<IN>){
    if(/^>/){
	my $gi =  (split(/\|/))[1];			#gi of fasta sequence
	say OUT "$gi\t$gihash{$gi}";
    }
}

#lucky draw
#my %gihash3;
#foreach (keys %gihash2) { 
#    my @array = @{$gihash2{$_}};
#    my $chosengi = $array[int(rand @array)];
#	$gihash3{$chosengi}++
#}
#say 'Theres ', scalar keys %gihash3, "chosen";

#my $in  = Bio::SeqIO->new(-file => "zcat out/sim.0106.input.gz |", -format => 'Fasta');
#my $out= IO::Zlib->new("out/sim.0107.out.gz", "wb9");
#
#while (my $seq = $in->next_seq ){
#    my $seqid 		= $seq->display_id;
#    my ($gi) =(split(/\|/, $seqid))[1];
##	key:giid	   value:refseq
#    if(exists $gihash3{$gi}){
#	say $out ">",$seqid,"|",$gihash{$gi};
#	say $out $seq->seq;
#    }
#}
#$out->close;
