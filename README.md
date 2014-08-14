#Part0: Init

##sim.0100 
Generate a ranked GENUS-ABUNDANDANCE list 

***output: namelist***
```
    taxon   total   rank (abundance)
    Caldilinea      25414.868178775 1
    Nitrospira      20661.7492551171        2
    Sorangium       18241.7170368186        3
    Mycobacterium   14892.9966161738        4
    Candidatus Accumulibacter       11906.5856237162        5
```

##sim.0010
 Choose 100 genomes based on refseq genome availablility
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

#Part2:	Remove related sequences from refseq

##sim.0020
 Determine the exact taxid of the genera 
 WHY? Redundancy in the NAME <-> TAXID 

##sim.0101 [INPUT:out/sim.0020.out.txt]

Takes chosen taxa list [data/topChosenGenera] and returns all taxa falling under the genera 
NOTE: searching by name returns taxa with same name but not under bacteria

***OUTPUT***
```
  originID     originName targetID                     targetName    rank
  1   104175 Oscillochloris   543045   uncultured Oscillochloris sp species
  2   104175 Oscillochloris   543045   uncultured Oscillochloris sp species
  3   104175 Oscillochloris   108012          Oscillochloris sp A19 species
  4   104175 Oscillochloris   108011           Oscillochloris sp BM species
  5   104175 Oscillochloris   104176      Oscillochloris trichoides species
  6   104175 Oscillochloris   765420 Oscillochloris trichoides DG-6 no rank
```

##sim.0103-105
   1. **0103** Takes the taxa gi list (from ncbi taxonomy), stores the targetID (TAXON ID) and stores the taxaID:gi hash table  [INPUT: sim.0101.out.txt]
   2. **0104** loops through refseq fasta files spits out 2 outputs [trimmed removed of GI] + [the removed sequences]
   3. **0105** batch job 104 across each individual gz file

##sim.0106
combines all the trimmed gz into 1

#Part2: Extract the genomes 

##sim.0200 [input sim.0010.out.txt]
	takes in refseq ID 
scans through refseq to extract the genome

##sim.0201
Replicates the abundances in the fasta file

##sim.0300
pushes the genomes into xc's scripts

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
