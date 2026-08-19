# RAXML

Maximum-likelihood phylogenetic inference outputs.

The best-scoring tree used by the dashboard, `COI_subsampled20_per_country_year.aligned.raxml.bestTree`,
now lives in [`../Dashboard/data/`](../Dashboard/data/) rather than here.

The app reads it at startup, and a Shiny deployment bundle only contains the app
directory — anything referenced with `../` would be missing once published. Keeping
the tree inside `Dashboard/data/` makes the app directory self-contained.

Re-running RAxML? Write the new `.bestTree` to `Dashboard/data/` under the same
filename, or update `TREE_PATH` at the top of `Dashboard/app.R`.
