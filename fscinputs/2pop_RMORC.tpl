//Number of population samples (demes)
2
//Population effective sizes (number of genes)
N1BOT$
N2BOT$
//Samples sizes and samples age
12
12
//Growth rates: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 0
0 0
//Migration matrix 1
0 M21$
M12$ 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
4 historical event
TBOT$ 0 0 0 RESIZE1$ 0 0
TBOT$ 1 1 0 RESIZE2$ 0 0
TISO$ 0 0 0 1 0 1
TDIV$ 1 0 1 RESIZEA$ 0 0
//Number sbof independent loci [chromosome]
1000000 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
DNA 100 0 0.0000000045
