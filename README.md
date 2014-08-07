#Functions of each script

##sim.0100
 From the relative abundance data choose the top100 genera [output: namelist]

##sim.0101
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

##sim.0107 (May be redundant; could use 103's hashtable)
  Outputs taxid for the refseq GIs belonging to the top100 


##sim.0109 
  Takes in the list of GIs with sequences and chooses one genome for each genera
