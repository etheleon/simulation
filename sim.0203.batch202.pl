#!/usr/bin/env perl 

use 5.012; 
use autodie;
use lib "/export2/home/uesu/perl5/lib/perl5";

#waiting##########################################
my $jobName = qx(qstat -u "uesu" | grep job6 | wc -l);
my $time =0;
while($jobName != 0) {  
$jobName = qx(qstat -u "uesu" | grep job6 | wc -l); 
say 'Waiting for the rest to finish. Time Elapsed: ',$time;
$time += 500;
sleep 500;
}
##################################################

my @refseq = </export2/home/uesu/db/refseq/arch2/*genomic.fna.gz>;
open my $output, ">", "out/sim.0203.out";
foreach my $fna (@refseq)
		{
#trim the fna.gz
	$fna =~ s/\.gz$//;
say $output "script/sim.0202.trimDB.pl $fna"
		}

my $numberofjobs = scalar @refseq;
system("qsub -b y -cwd -V -t 1-$numberofjobs -tc 60 /export2/home/uesu/local/bin/rasengan out/sim.0203.out");
