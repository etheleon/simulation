#!/usr/bin/env perl 

use 5.012;
use lib "/export2/home/uesu/perl5/lib/perl5"; 
use autodie;
use Bio::SeqIO;
use IO::Zlib;
use Storable;

die "USAGE: $0 <genomes.fna> \n" unless $#ARGV == 0;

#Stores the gi2taxid
my %gihash = %{retrieve("out/sim.0103.out.pdo")};

my $in  = Bio::SeqIO->new(-file => "$ARGV[0]", -format => 'Fasta');
my $out = Bio::SeqIO->new(-file => ">out/sim.0106.out.fna", -format => 'Fasta');
say "#creating gi-taxa list";

while (my $seq = $in->next_seq ){
    my $seqid 		= $seq->display_id;
    my ($gi) =(split(/\|/, $seqid))[3];
#	key:giid	   value:refseq
    if(exists $gihash{$gi}){
    	$out -> write_seq($seq);
    }
}
$out->close;
