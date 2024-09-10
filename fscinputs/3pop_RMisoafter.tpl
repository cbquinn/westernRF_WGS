//Number of population samples (demes)
3
//Population effective sizes (number of genes)
N1BOT$
N2BOT$
N3BOT$
//Samples sizes and samples age
12
12
10
//Growth rates: negative growth implies population expansion
0
0
0
//Number of migration matrices : 0 implies no migration between demes 
4
//Migration matrix 0: no migration (initial state) 
0 0 0
0 0 0
0 0 0
//Migration matrix 1: migration between pop 1 and pop 2 only
0 0 0
0 0 Morclas$
0 Mlasorc$ 0
//Migration matrix 2: migration between pop 0 and pop 1 and op 0 and pop 2
0 Mrmorc$ Mrmlas$
Morcrm$ 0 0
Mlasrm$ 0 0
//Migration matrix 2: migration between pop 0 and pop 1 only
0 Mrmsn$ 0
Msnrm$ 0 0
0 0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
8 historical event
30 0 0 0 RESIZE1$ 0 0
30 1 1 0 RESIZE2$ 0 0
30 2 2 0 RESIZE3$ 0 0
TISO.las$ 0 0 0 1 0 1
TISO.rm$ 0 0 0 1 0 2
TDIV.orclas$ 2 1 1 RESIZEB$ 0 3
TDIV.orclas$ 0 0 0 RESIZEC$ 0 3
TDIV.rmsn$ 1 0 1 RESIZEA$ 0 0
//Number sbof independent loci [chromosome]
1000000 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
DNA 100 0 0.0000000045
