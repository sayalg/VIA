# Dashboard/data — generated, do not edit

Everything here is written by `Rscript Tools/build_app_cache.R` (run from the
repository root). To change what the dashboard shows, edit `Input_Files/` and
rebuild — see "Updating the data" in the top-level README.

| File | Contents |
|---|---|
| `genes.rds` | One entry per folder in `Input_Files/genes/`: the map table, haplotype palette, year range, StrainHub network, country connectivity network (with precomputed layout), and tree-figure path. |
| `map_layers.rds` | Country native/invasive status polygons, and the railway and shipping-route geometry pre-serialised to GeoJSON. |
| `pubmed_snapshot.rds` | Seed for the Literature tab, so the first visitor isn't left waiting 20–40 s for PubMed. If the build can't reach PubMed it leaves the previous snapshot in place. |

The build also copies each gene's `tree.svg` to `Dashboard/www/genes/<GENE>/`.

These files are committed so that `Dashboard/` can be deployed straight from the
repository.
