#!/usr/bin/env perl

use strict;
use v5.10;
use lib "/export2/home/uesu/perl5/lib/perl5"; 
use autodie;
use Storable;

#Stores the GI to be removed ie. belonging to the top100 genera
my %taxahash;
open TOBEDELETED, "/export2/home/uesu/simulation_fr_the_beginning/out/sim.0101.out.txt";
while(<TOBEDELETED>){
unless(/^#/ || /^taxid/){
    chomp;
my ($taxa, $parent) =  (split(/\t/))[0,3];
$taxahash{$taxa}++ unless $parent =~ /33057/;	#Keep Thauera
}
}
close TOBEDELETED;
say "#Finished storing condemned taxa. Total",scalar keys %taxahash;

#assigns the gi belonging to condemned taxa to hash (key::taxid)
my %gihash;
open GITAX, "zcat /export2/home/uesu/db/taxonomy_latest/gi_taxid_nucl.dmp.gz | ";
while(<GITAX>){
    chomp;
    my($gi, $taxid) = split(/\t/);
#    say "$gi\t$taxid";
    $gihash{$gi}++ if (exists $taxahash{$taxid});
}
close GITAX;
say "#Finished stored condemned gi. Total",scalar keys %gihash;

store \%gihash, "out/sim.0103.out.pdo";

