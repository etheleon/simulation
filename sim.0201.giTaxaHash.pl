#!/usr/bin/env perl

use strict;
use v5.10;
use lib "/export2/home/uesu/perl5/lib/perl5"; 
use autodie;
use Storable;

die "$0 <out/sim.0200.out.txt> \n" unless $#ARGV == 0;

#Stores the GI to be removed ie. belonging to the top100 genera
my %taxahash;
#open TOBEDELETED, "/export2/home/uesu/simulation_fr_the_beginning/out/sim.0101.out.txt";
open TOBEDELETED, $ARGV[0];
while(<TOBEDELETED>){
unless(/^#/ || /^originID/){
    chomp;
    my ($parent, $taxa) =  (split(/\t/))[0,2];
    $taxahash{$taxa}=$parent unless $parent =~ /33057/;	#Keep Thauera
}
}
close TOBEDELETED;
say "#Finished storing condemned taxa. Total ", scalar keys %taxahash;

#assigns the gi belonging to condemned taxa to hash (key::taxid)
my %gihash;
open GITAX, "zcat /export2/home/uesu/db/taxonomy_latest/gi_taxid_nucl.dmp.gz | ";
while(<GITAX>){
    chomp;
    my($gi, $taxid) = split(/\t/);
#    say "$gi\t$taxid";
    $gihash{$gi} = $taxahash{$taxid} if (exists $taxahash{$taxid});	#links GI with parentGenus
}
close GITAX;
say "#Finished stored condemned gi. Total",scalar keys %gihash;
store \%gihash, "out/sim.0201.out.pdo";
