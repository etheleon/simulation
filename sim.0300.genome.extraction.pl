#!/usr/bin/env perl  

use lib "/export2/home/uesu/perl5/lib/perl5"; 
use 5.012;
use autodie;
use Bio::SeqIO;

die "$0 <sim.0010.out.txt>\n" unless $#ARGV == 0;

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

##################################################
#Part2: reads all the files individually
##################################################

my $out = Bio::SeqIO->new(-file => ">out/sim.0200.out.fna", -format => 'Fasta');
my @files = <"/export2/home/uesu/db/refseq/arch2/*">; 

foreach (@files){ 
    my $in  = Bio::SeqIO->new(-file => "zcat $_ |", -format => 'Fasta');
    while (my $seq = $in->next_seq ){
    	my $seqid = $seq->display_id;
    	my ($refid) =(split(/\|/, $seqid))[3];
	    if(exists $ref{$refid}) { 
	    	$seqid = 'genustaxaid|'.$ref{$refid}.'|'.$seqid;
	    	$seq->display_id($seqid);
	    	$out->write_seq($seq);
		delete $ref{$refid};
    	    }
    }
    $in->close;
}
$out->close;


#Which entries have been removed
open my $report, '>', 'out/sim.0200.out.report';
say $report "refseq\ttaxid";
foreach (keys %ref) { say $report "$_\t$ref{$_}"}
close $report;
