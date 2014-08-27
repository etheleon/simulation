#!/usr/bin/env perl  

use lib "/export2/home/uesu/perl5/lib/perl5"; 
use 5.012;
use autodie;
use Bio::SeqIO;

die "$0 <sim.0101.out.txt> <out/sim.0103.chosen.txt>\n" unless $#ARGV == 1;

##################################################
#Part1: reads in the NC of the chosen genomes
##################################################
my %ref;
open my $in, 	'<', 	$ARGV[0];

while(<$in>) {
    unless($. == 1){
    	chomp;
    	my ($taxid, $refseqid) = (split(/\t/))[1,3];
    	my @refseq = split /,/, $refseqid;
    	if( scalar @refseq > 1) { 
	    foreach(@refseq) {$ref{$_} = $taxid}
    	}else{$ref{$refseqid}=$taxid}
    }
}
close $in;
say scalar keys %ref, " refIDs of complete genomes taken in"; 

open my $incomplete, '<', $ARGV[1];
while(<$incomplete>) { 
unless ($. == 1) 
	{ 
	    	chomp;
	    	my ($taxid, $refseqid) = (split(/\t/))[1,3];
		$ref{$refseqid} = $taxid;	#
	}
}
close $incomplete;
say scalar keys %ref, " refIDs of complete and incomplete genomes taken in";

##################################################
#Part2: reads all the files individually
##################################################

my $out = Bio::SeqIO->new(-file => ">out/sim.0300.out.fna", -format => 'Fasta');
my @files = <"/export2/home/uesu/db/refseq/arch2/*genomic.fna.gz">; 

open my $outt, '>', "out/sim.0300.out.length.txt";
say $outt join "\t", qw/refseqID genus length/;#header 
foreach (@files){ 
    my $in  = Bio::SeqIO->new(-file => "zcat $_ |", -format => 'Fasta');
    while (my $seq = $in->next_seq ){
    	my $seqid = $seq->display_id;
    	my ($refid) =(split(/\|/, $seqid))[3];
	    if(exists $ref{$refid}) { 
	    	$seqid = 'genustaxaid|'.$ref{$refid}.'|'.$seqid;
	    	$seq->display_id($seqid);
	    	$out->write_seq($seq);
		my $length = $seq->length;
		say $outt join "\t", $refid, $ref{$refid}, $length; 
		delete $ref{$refid};
    	    }
    }
    $in->close;
}
$out->close;


#Which entries have not been removed
open my $report, '>', 'out/sim.0300.out.report';
say $report "refseq\ttaxid";
foreach (keys %ref) { say $report "$_\t$ref{$_}"}
close $report;
