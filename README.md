#Part1: Init

***sim.0100.genometaxid.r***
---

Function: Generate a ranked GENUS-ABUNDANDANCE list 

|Input | Output |
|------|--------|
|abundance table| data/top500genera.gDNA|

```
    taxon   total   rank (abundance)
    Caldilinea      25414.868178775 1
    Nitrospira      20661.7492551171        2
    Sorangium       18241.7170368186        3
    Mycobacterium   14892.9966161738        4
    Candidatus Accumulibacter       11906.5856237162        5
```

***sim.0101.chooseGenomes.r***


THINGS TO DO: 
1. need to include up to 500 genera
2. Those with complete genomes take complete genomes; for those without, 
  1. identify all of the children 
  2. comb through REFSEQ for GIs associated with the child taxids 
  3. and summarise 
  combined length of the sequences with each of the taxas  
  Somehow need to figure out how to choose the genomes
	
 Choose for genera with complete genomes, based on refseq genome availablility (refer to dl-ed table)
 Gives the following outputs:
	1. the selected genera (name) to be removed [topChosenGenera.txt]
```
    Streptomyces    9883.98971459188        6
    Anaerolinea     8975.14623219175        7
    Burkholderia    8515.90087908299        8
    Thauera 6212.32458541424        9
    Plesiocystis    6043.42026585336        10
```
	2. the refseq IDs of the randomly chosen genomes [sim.0010.out.txt]

NOTE (not enough that we include only genera with complete genomes)
THINGS TO DO need script to find the missing ones and 

***sim.0102.combThruGIs.pl***
[input sim.0101.missing.txt]
1. Takes in taxids of genera without complete genomes & finds all children taxa
2. Builds hash of childrentaxa::gi
3. Loops through refseq and takes the sequence refseqID 

OUTPUT
```
         gi  taxid parentGenus             refid
         1 254971245 263906      499551       NR_028163.1
         2 631251601 547055      499551       NR_112799.1
         3  67926185 165597      263510 NZ_AADV02000295.1
         4  67926203 165597      263510 NZ_AADV02000300.1
         5  67926213 165597      263510 NZ_AADV02000303.1
         6  67926222 165597      263510 NZ_AADV02000306.1
```
gives output for the total number of sequences per taxa


***sim.0103.summary.r***
summarises the 
	# of sequences / parent genera
	# of sequences / taxa
	




#Part2:	Remove related sequences from refseq

***sim.0200.get.the.species.underneath*** 
[INPUT:out/sim.0020.out.txt]
Takes chosen taxa list [data/topChosenGenera] and returns all taxa falling under the genera 
NOTE: searching by name returns taxa with same name but not under bacteria

[OUTPUT]
```
  originID     originName targetID                     targetName    rank
  1   104175 Oscillochloris   543045   uncultured Oscillochloris sp species
  2   104175 Oscillochloris   543045   uncultured Oscillochloris sp species
  3   104175 Oscillochloris   108012          Oscillochloris sp A19 species
  4   104175 Oscillochloris   108011           Oscillochloris sp BM species
  5   104175 Oscillochloris   104176      Oscillochloris trichoides species
  6   104175 Oscillochloris   765420 Oscillochloris trichoides DG-6 no rank
```

***sim.0201-203***
   1. **0201.giTaxaHash.pl** Takes the taxa gi list (from ncbi taxonomy), stores the targetID (TAXON ID) and stores the taxaID:gi hash table  [INPUT: sim.0101.out.txt]
   2. **0202.trimDB** loops through refseq fasta files spits out 2 outputs [trimmed removed of GI] + [the removed sequences]
   3. **0203** batch job 202 across each individual gz file

***sim.0204.combinedTrimmed.sh***
combines all the trimmed gz into 1

#Part3: Extract the genomes 

***sim.0300***
[input sim.0010.out.txt]
	takes in refseq ID 
scans through refseq to extract the genome

#Part4: read creation
***sim.0300***
Pushes the genomes into xc's scripts

Counting the number of sequences 
```bash
   cd /export2/ulu_pandan2011/
   for i in `ls`; do cat $i | echo $i\t$((`wc -l`/4)); done;

s_1_1.filtered.fastqt 172976284
s_1_2.filtered.fastqt 172976284

s_2_1.filtered.fastqt 164554783
s_2_2.filtered.fastqt 164554783

s_3_1.filtered.fastqt 163325133
s_3_2.filtered.fastqt 163325133

s_4_1.filtered.fastqt 166683965
s_4_2.filtered.fastqt 166683965

s_5_1.filtered.fastqt 166518064
s_5_2.filtered.fastqt 166518064

s_6_1.filtered.fastqt 164701017
s_6_2.filtered.fastqt 164701017

s_7_1.filtered.fastqt 170304759
s_7_2.filtered.fastqt 170304759

s_8_1.filtered.fastqt 182780732
s_8_2.filtered.fastqt 182780732
```
