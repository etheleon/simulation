1: From the relative abundance data choose the top100 genera [output: namelist]

2: sim.0101 Takes name list and returns all taxa related to the top100 list
  originID     originName targetID                     targetName    rank
  1   104175 Oscillochloris   543045   uncultured Oscillochloris sp species
  2   104175 Oscillochloris   543045   uncultured Oscillochloris sp species
  3   104175 Oscillochloris   108012          Oscillochloris sp A19 species
  4   104175 Oscillochloris   108011           Oscillochloris sp BM species
  5   104175 Oscillochloris   104176      Oscillochloris trichoides species
  6   104175 Oscillochloris   765420 Oscillochloris trichoides DG-6 no rank

3: sim.0103 
   Takes the list, stores the targetID and stores a hash table [key=gi]

4: sim.0104+sim.0105 
   loops trhough refseq fasta files spits out 2 outputs [trimmed] + [removed]
   trimmed is the file removed of the gi, and removed are the removed sequences

5. sim.0106 
   Combines the removed gi into 1 file. 

6. sim.0107
   Outputs the files in the combined removed list gi, taxa to file sim.0107.out.txt

7  sim.0108
   Output
