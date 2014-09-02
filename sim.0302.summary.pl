#!/usr/bin/env perl

use warnings; 
use strict;
use v5.10;

use lib "/export2/home/uesu/perl5/lib/perl5";
use autodie;
use Bio::SeqIO;

say "ID\tLENGTH";
my  $seqIN= Bio::SeqIO->new( -format => 'Fasta', -file => 'out/sim.0300.combined.fna');
	    while (my $seq = $seqIN->next_seq ){
		my $id  	= $seq->display_id;
		my $length 	= $seq->length;
		say "$id\t$length";
	    }
$seqIN->close;
