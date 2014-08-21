#!/usr/bin/env perl 

use 5.012; 
use autodie;
use lib "/export2/home/uesu/perl5/lib/perl5";

my @refseq = </export2/home/uesu/db/refseq/arch2/*>;
open my $output, ">", "out/sim.0105.out";
foreach my $fna (@refseq)
		{
		say $output "qsub -V -cwd -b y script/sim.0104.builddb.pl $fna"
		}

my $numberofjobs = scalar @refseq;
system("qsub -b y -cwd -V -t 1-$numberofjobs -tc 60 /export2/home/uesu/local/bin/rasengan out/sim.0105.out");