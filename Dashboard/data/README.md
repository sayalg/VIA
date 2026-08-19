# Dashboard/data

Runtime data for `app.R`. Everything the app reads at startup lives here so that
the `Dashboard/` directory is fully self-contained and can be published as-is.

**Do not move these files out of this directory.** A Shiny deployment bundle
(`rsconnect::deployApp("Dashboard")`, or a Shiny Server app directory) contains
only the app directory. A path like `../Input_Files/...` resolves fine locally
but is missing once published.

| File | Size | Read by | Produced by |
|---|---|---|---|
| `FinalCOI_metadata.csv` | 163 KB | `SAMPLES_PATH` — sample table, location cleaning | Curated by hand |
| `FinalCOI_metadata_w_coor.csv` | 171 KB | `METADATA_PATH` — map markers + both networks | `fill_missing_coords()` in `app.R`, run once offline |
| `MTM_INVASIVE_VECTOR_SPECIES_20251205.xlsx` | 189 KB | `INVASIVE_PATH` — country native/invasive layer (sheet `Data`, columns `COUNTRY_NAME`, `INVASIVE_STATUS`) | Malaria Threat Map export |
| `COI_subsampled20_per_country_year.aligned.raxml.bestTree` | 15 KB | `TREE_PATH` — StrainHub transmission network | RAxML (see `../../RAXML/`) |
| `africa_railways.rds` | 5.6 MB | `RAILWAYS_PATH` — railway overlay (optional) | `../../Railway_Extraction/prepare_africa_railways.R` |
| `shipping_routes.rds` | 7.1 MB | `SHIPPING_PATH` — shipping-route overlay (optional) | `../../Shipping_Extraction/prepare_shipping_routes.R` |

The two `.rds` files are optional: if either is absent the app starts with a
warning and hides the corresponding map layer. The other four are required and
`app.R` stops at startup with an explicit message if any is missing.

Raw sequence data (`.fasta`) is **not** here — it is analysis input, not app
input, and stays in [`../../Input_Files/`](../../Input_Files/).
