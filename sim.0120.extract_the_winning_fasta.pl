#!/usr/bin/env perl 

use v5.10; 
use strict;

use lib "/export2/home/uesu/perl5/lib/perl5";
use Bio::SeqIO;
use autodie;

my %gihash;
open IN, "out/sim.0109.out.txt"; 
while(<IN>) {
    unless ($. == 1) { 
	chomp;
	my $gi = (split(/\s+/))[0];
	$gihash{$gi}++;
    }
}

my $in 	= Bio::SeqIO->new(-file => "zcat out/sim.0106.input.gz | ", -format => 'Fasta');
my $out = Bio::SeqIO->new(-file => ">out/sim.0120.output.fna ", -format => 'Fasta');

while(my $seq = $in -> next_seq) { 
    my $seqid 		= $seq->display_id;
    my ($gi) =(split(/\|/, $seqid))[1];
#	key:giid	   value:refseq
	if(exists $gihash{$gi}) { 
	    $out -> write_seq($seq);
	}
}
