#!/usr/bin/env perl 

use strict;
use v5.10;
use Storable;
use autodie;

say "Loading the full gi-taxa hash";
my %gihash = %{retrieve("out/sim.0103.out.pdo")};

my %gihash2;
say "Loading 2nd gi-taxa hash";
open INPUT, "out/sim.0107.out.txt";
while(<INPUT>){
    chomp;
    my ($gi, $taxid) = split(/\t/);
    $gihash2{$gi} = $taxid;
}

say "Both loaded";

my @this_not_that = ();
foreach (keys %gihash) {
        say "$_\tnoting2" unless exists $gihash2{$_};
}

my @this_not_that2 = ();
foreach (keys %gihash2) {
        say  "$_\tnotin" unless exists $gihash{$_};
}
