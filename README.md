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
---

|Input | Output | Comments | 
|------|--------|----------|
|data/genomereports/prokaryotes.txt|| downloaded from ftp://ftp.ncbi.nih.gov/genomes/GENOME_REPORTS updated daily|
|data/top500genera.gDNA|||
||out/sim.0101.out.txt||
||out/sim.0101.out2.txt|genera to be included (8 of which do not have any refseq sequences)|
||data/topChosenGenera.txt||
||out/sim.0101.missing.txt||

Filters:
1. refseqID 
2. gapless chromosome

sim.0101.out.txt::265 genomes from 265 genera 
```
origin  taxid   taxon   Chromosomes.RefSeq
1284663 1578    Lactobacillus   NC_020229.1
231023  286     Pseudomonas     NC_017986.1
680198  1883    Streptomyces    NC_013929.1
1234365 209     Helicobacter    NC_018937.1
1032845 780     Rickettsia      NC_015866.1
488222  1301    Streptococcus   NC_012466.1
1005048 202907  Collimonas      NC_015856.1
```
data/topChosenGenera.txt::name/total.G.Abundance/rank
```
    taxon           total             	 rank
    Streptomyces    9883.98971459188        6
    Anaerolinea     8975.14623219175        7
    Burkholderia    8515.90087908299        8
    Thauera 6212.32458541424        9
    Plesiocystis    6043.42026585336        10
```
out/sim.0101.missing.txt
```
taxid
1234
233189
191767
162027
2349
354354
113
```

Summary
1. 	Those with complete genomes take complete genomes; for those without (see next script)
  * Choose for genera with complete genomes, based on refseq genome availablility (ftp://ftp.ncbi.nih.gov/genomes/GENOME_REPORTS/)

***sim.0102.combThruGIs.pl***
[input sim.0101.missing.txt]

|Input | Output | Comments | 
|------|--------|----------|
|out/sim.0101.missing.txt	|			|					|
|				|out/sim.0102.out	| gi-child-genus-refseqID		|
|				|out/sim.0102.out2	| genus - bplength 			|
|				|out/sim.102.output.fna	| FASTA sequences of children taxa	|

out/sim.0102.out
```
gi      taxid   parentGenus     refid
254971245       263906  499551  NR_028163.1
631251601       547055  499551  NR_112799.1
67926185        165597  263510  NZ_AADV02000295.1
67926203        165597  263510  NZ_AADV02000300.1
```

  1. Identify children of genera without complete genomes
  2. Comb through REFSEQ for GIs associated with the child taxids 


***sim.0103.summary.r***
Some text:

|Input | Output | Comments | 
|------|--------|----------|
|out/sim.0102.out | 		|				|
|out/sim.0102.out2|		|				|
||out/sim.0103.chosen.txt	|refseqIDs of chosen taxa-genus	|
||out/sim.0103.summary1.pdf	|				|

summarises the 
	# of sequences / parent genera
	# of sequences / taxa
	
#Part2:	Remove related sequences from refseq

***sim.0200.get.the.species.underneath.pl*** 

sometext 

|Input | Output | Comments | 
|------|--------|----------|
|out/sim.0101.out2.txt||list of genera - taxid (with or without genomes)|


Takes chosen taxa list and returns all taxa falling under the genera (genera with or without complete genomes)
NOTE: searching by name returns taxa with same name but not under bacteria

out/sim.0200.out.txt
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
   1. **0201.giTaxaHash.pl** Takes the taxa gi list (from ncbi taxonomy), stores the targetID (TAXON ID) and stores the taxaID:gi hash table  [INPUT: sim.0200.out.txt]
   2. **0202.trimDB** loops through refseq fasta files spits out 2 outputs [trimmed removed of GI] + [the removed sequences] [input:sim]
   3. **0203** batch job 202 across each individual gz file

#Part3: Extract the genomes 
```
tail -n +2 out/sim.0101.out.txt
tail -n +2 out/sim.0103.chosen.txt | perl -aln -F"\t" -E 'say $F[1]' | sort | uniq| wc -l
251+162
```

***sim.0300***
	takes in refseq ID of both complete + incomplete
	scans through refseq to 
1. extract the genome
2. refseq seq of genera without genomes and concatenate into1


#Part4: read creation
***sim.0400***
needs out/sim.0101.out2.txt (abundance)

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
