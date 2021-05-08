# phylogeographySARS-CoV-2
parsimony phylogeography analyses for SARS-CoV-2

©EMBL-European Bioinformatics Institues, 2021

This is the code for performing phylogeography analyses of SARS-CoV-2 VOCs data to infer introductions into the UK and clade sizes. The code is highly ad-hoc for the considered project and hard-coded, but the principles in the code could be more generally useful, let me know if you are interested.

The idea behind the approach used here is to consider a fixed phylogeny, inferred by maximum likelihood with FastTree2, and then estimate introductions into the UK using parsimony, but accounting for all possible splits of multifurcation nodes in the tree, and acounting for all possible parsimonious migration histories. The approach will calculate minimum and maximum number of introductions within parsimonious migration histories and give an example of clade sizes for one such parsimonious migration history. The approach for parsimony introduction inference is very fast as it uses dynamic programming along the phylogenetic tree, so the cost is linear in tree size (that is, linear in the number of samples considered).
