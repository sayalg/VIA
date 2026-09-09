# Dashboard/data

Runtime data for `app.R`. Everything the app reads at startup lives here so that
the `Dashboard/` directory is fully self-contained and can be published as-is.

**Do not move these files out of this directory.** A Shiny deployment bundle
(`rsconnect::deployApp("Dashboard")`, or a Shiny Server app directory) contains
only the app directory. A path like `../Input_Files/...` resolves fine locally
but is missing once published.

## Contents

| File | Size | Read by | Produced by |
|---|---|---|---|
| `FinalCOI_metadata.csv` | 163 KB | `SAMPLES_PATH` — sample table, location cleaning | Curated by hand |
| `FinalCOI_metadata_w_coor.csv` | 171 KB | `METADATA_PATH` — map markers | `fill_missing_coords()` in `app.R`, run once offline |
| `MTM_INVASIVE_VECTOR_SPECIES_20251205.xlsx` | 189 KB | `INVASIVE_PATH` — country native/invasive layer (sheet `Data`, columns `COUNTRY_NAME`, `INVASIVE_STATUS`) | Malaria Threat Map export |
| `COI_subsampled20_per_country_year.aligned.raxml.bestTree` | 15 KB | `TREE_PATH` — fallback network build only | RAxML (see `../../RAXML/`) |
| `africa_railways.rds` | 176 KB | `RAILWAYS_PATH` — railway overlay | **`Tools/build_app_cache.R`** |
| `shipping_routes.rds` | 58 KB | `SHIPPING_PATH` — shipping overlay | **`Tools/build_app_cache.R`** |
| `networks.rds` | 1 KB | Both visNetwork panels | **`Tools/build_app_cache.R`** |
| `pubmed_snapshot.rds` | 233 KB | Seeds the Literature tab | **`Tools/build_app_cache.R`** |

## Generated caches

The last four are **build artefacts, not source data**. Regenerate them from the
repository root whenever the underlying data changes:

```bash
Rscript Tools/build_app_cache.R
```

They exist because the app was shipping ~75 MB of map geometry to every browser,
which is what made it slow to load and prone to crashing on a small instance:

- **`shipping_routes.rds`** — the raw extract repeats each from–to pair once per
  year it was active: 34,735 features where only 1,774 are distinct, all drawn on
  top of one another. The cache stores each geometry once plus a
  `route_id → year` index. 57 MB of GeoJSON becomes 2.6 MB, with no visible
  difference. Deliberately *not* simplified — searoute paths hug coastlines, and
  shaving vertices can walk a sea lane across land.
- **`africa_railways.rds`** — three things happen here, and together they take the
  layer from ~11 s of server time per session to ~2.6 s:
  1. **Short segments are dropped** (`RAILWAY_MIN_SEGMENT_M`, default 500 m).
     OpenStreetMap splits a railway at every attribute change, so the median raw
     segment is 355 m and a quarter are under 100 m — sidings and yard tracks,
     invisible at the zooms this map is used at. They dominate the cost, because
     every segment keeps at least two vertices no matter how hard you simplify.
     At 500 m the layer keeps **94.6% of the network by length** for 40% of the
     vertices.
  2. **`st_combine`, not `st_union`.** `summarise()` on an sf object defaults to
     `st_union`, which nodes the lines at every intersection — that alone was
     adding ~50,000 vertices for no visible benefit.
  3. Simplified to 500 m and merged to one MULTILINESTRING per railway type.

  Net: 120,801 → 32,920 vertices, 2.7 MB → 793 KB on the wire. The trade is that
  individual railway names are gone from the hover label; the type is still
  shown, and the legend is unchanged.
- **`networks.rds`** — both visNetwork node/edge tables, so no session pays for
  the parsimony reconstruction or the centrality maths.
- **`pubmed_snapshot.rds`** — seed results so the first visitor after a process
  start isn't left waiting 20–40 s for the live NCBI fetch.

If a cache is missing the app still runs: the two network panels fall back to
computing at runtime, and the two map layers hide themselves with a startup
warning.

## Raw inputs

The raw pipeline artefacts these are built from live outside the app directory,
in [`../../Input_Files/`](../../Input_Files/) — `africa_railways_raw.rds` and
`shipping_routes_raw.rds`, along with the `.fasta` sequence data, which is
analysis input rather than app input.
