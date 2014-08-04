#!/usr/bin/env perl 

use strict;
use lib "/export2/home/uesu/perl5/lib/perl5"; 
use v5.10;
use autodie;
use Bio::SeqIO;
use IO::Zlib;
use Storable;

die "USAGE: $0 <refseq> leave the file extension out\n" unless $#ARGV == 0;

my %gihash = %{retrieve("out/sim.0103.out.pdo")};

##################################################
#Summary theres a total of 1260572 sequences
##################################################

#NOTE if fasta has gi belonging to the condemned list dont write to db but print to stdout

my $in  = Bio::SeqIO->new(-file => "zcat $ARGV[0].gz |", -format => 'Fasta');

my $outputfilename = "$ARGV[0]"."_trimmed.gz";
my $outputfilename2 = "$ARGV[0]"."_removed.gz";
my $out= IO::Zlib->new("$outputfilename", "wb9");
my $out2= IO::Zlib->new("$outputfilename2", "wb9");
#
say "#Writing new DB and remove condemned";
#
while (my $seq = $in->next_seq ){
    my $seqid 		= $seq->display_id;
    my ($gi, $refid) =(split(/\|/, $seqid))[1,3];
#	key:giid	   value:refseq
    if(!exists $gihash{$gi}){
	say $out ">",$seqid;
	say $out $seq->seq;
    }else{
	say "foundone";
 	say $out2 ">",$seqid;
	say $out2 $seq->seq;
    }
}
$out->close;
