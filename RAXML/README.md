# RAXML

Maximum-likelihood phylogenetic inference outputs.

The best-scoring tree the dashboard uses lives with the rest of its gene's
inputs, at `Input_Files/genes/<GENE>/tree.newick` (for COI, formerly
`COI_subsampled20_per_country_year.aligned.raxml.bestTree`).

Re-ran RAxML? Copy the new `.bestTree` over that `tree.newick`, re-export
`tree.svg` from iTOL if the figure should change, and run
`Rscript Tools/build_app_cache.R` from the repository root.
