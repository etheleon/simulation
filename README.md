#Functions of each script

##sim.0100 
 (simple)
 From the relative abundance data choose the top100 genera [output: namelist]

##sim.0101 
(submission requires NEO$J)
Takes name list and returns all taxa related to the top100 list
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
   1. **0103** Takes the taxa gi list, stores the targetID (TAXON ID) and stores the taxaID:gi hash table 
   2. **0104** loops through refseq fasta files spits out 2 outputs [trimmed removed of GI] + [the removed sequences]
   3. **0105** batch job 104 across each individual gz file

##sim.0106 
   Combines the removed gi into 1 file. 

##sim.0107 
  May be redundant; could use 103's hashtable, using checking.pl
  Outputs taxid for the refseq GIs belonging to the top100 

##sim.0109 
  Takes in the list of GIs with sequences and chooses one genome for each genera


##sim.0120
  Takes the list of selected GIs and extracts the sequences

##sim.0120
  I need to know how the refseq sequences are taxa

##sim.0121
   Cal the abundance of the genera and corresponding the number of sequences



   ```
   cd /export2/ulu_pandan2011/
   for i in `ls`; do cat $i | echo $i\t$((`wc -l`/4)); done;

s_1_1.filtered.fastqt172976284
s_1_2.filtered.fastqt172976284
s_2_1.filtered.fastqt164554783
s_2_2.filtered.fastqt164554783
s_3_1.filtered.fastqt163325133
s_3_2.filtered.fastqt163325133
s_4_1.filtered.fastqt166683965
s_4_2.filtered.fastqt166683965
s_5_1.filtered.fastqt166518064
s_5_2.filtered.fastqt166518064
s_6_1.filtered.fastqt164701017
s_6_2.filtered.fastqt164701017
s_7_1.filtered.fastqt170304759
s_7_2.filtered.fastqt170304759
s_8_1.filtered.fastqt182780732
s_8_2.filtered.fastqt182780732
```
