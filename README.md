# VIA: Vector Invasion Atlas

An interactive R Shiny dashboard showing *Anopheles stephensi* haplotype distribution on a map (pie chart per collection site), phylogenetic trees, phylogeographic dispersal and haplotype-sharing networks, and a live hub of *An. stephensi* literature.

---

## Updating the data, or adding a gene

Everything the dashboard shows is rebuilt by **one command**, run from the repository root:

```bash
Rscript Tools/build_app_cache.R
```

It takes about 20 seconds. Afterwards `Dashboard/` is ready to run or deploy.

### Add new samples to an existing gene

1. Add rows to `Input_Files/genes/<GENE>/metadata.csv`.
   Leave `Latitude`/`Longitude` blank if you don't have them — the build looks them up from `City.Town`, `County.State` and `Country` via OpenStreetMap.
2. If the tree changed, replace `tree.newick` (and `tree.svg`, re-exported from iTOL) in the same folder.
3. Run the build.

### Add a new gene

1. Export the gene's sample sheet from Excel as CSV, saved as `metadata.csv`. Use the same columns as the other genes' sheets, plus the gene's own result column named **`<Gene> Genotype`** (or `<Gene> Haplotype`) — e.g. `P47Rec Genotype`.
2. Put it in a folder named after the gene: `Input_Files/genes/P47rec/metadata.csv`. **The folder name must match the start of that column** (capitals don't matter), because that's how the build knows which column to colour the map by. A `COI Haplotype` column, or any other extra column, is carried along and doesn't interfere.
3. Optionally add, in the same folder:

   | File | What it drives |
   |---|---|
   | `tree.newick` | StrainHub phylogeographic dispersal network (RAxML best tree; tip labels must match `FASTA ID`) |
   | `tree.svg` | Phylogenetic Tree tab (export it from iTOL) |

   Without them, those two panels show a note saying what's missing.
4. Run the build. The gene appears in the dashboard's gene dropdown, named after its folder, with the map legend titled after its column.

The sheet needs at least `Year`, `Country`, `Latitude`, `Longitude`, `City/Town`, `County/State` and the genotype column, plus `FASTA ID` if you add a tree. Excel-style headers (`City/Town`, `FASTA ID`) and R-style headers (`City.Town`, `FASTA.ID`) both work.

### What the build does for you

- **Missing coordinates** — first copied from any other gene's sheet that has the same `Sample ID` (so a sample sits in the same place under every gene), then looked up once per distinct address on OpenStreetMap and saved to `Input_Files/geocode_cache.csv`. Commit that file; reruns then make no network calls.
- **The same sample listed twice** (same Sample ID, genotype, country, town and year) — counted once, keeping the more complete row.
- **Misaligned FASTA IDs** — FASTA IDs end in their country (`Kassala8_Sudan`); if that disagrees with the row's Country, the build warns that the column has slipped, and refuses to build a StrainHub network from it.
- **Missing columns, an ambiguous genotype column, or other duplicates** — stops or warns, naming the gene and what to fix.
- **Samples it can't place** — warns and lists the addresses, so you can fix the spelling or type the coordinates in.
- Rows located at `Laboratory` are kept for the networks but left off the map.

### Not covered by the build

These produce the build's inputs, need their own tools, and only need rerunning when their sources change:

- **Alignment and variant calling** — `Scripts/*.sh` (MAFFT etc.).
- **RAxML** to produce `tree.newick`, and **iTOL** to draw `tree.svg`.
- **Railway and shipping extracts** — `Railway_Extraction/` and `Shipping_Extraction/` produce `Input_Files/*_raw.rds` (slow: Geofabrik downloads, the 40 MB UNCTAD table, Python `searoute`).

---

## Required R Packages

To **run** the dashboard:

```r
install.packages(c(
  "shiny", "bslib", "shinycssloaders", "leaflet", "leaflet.minicharts",
  "dplyr", "tidyr", "sf", "visNetwork", "stringr", "httr", "xml2"
))
```

To **rebuild the data** you also need the following; the build script checks for them and prints the exact `install.packages()` line for anything missing:

```r
install.packages(c(
  "readxl", "readr", "rnaturalearth", "rnaturalearthdata", "jsonlite", "tidygeocoder",
  "igraph", "ape", "tidytree", "castor", "hash", "plyr", "tibble", "magrittr"
))
```

Everything is on CRAN. Do **not** install the `strainhub` package: its `DESCRIPTION` declares none of its 27 imports, so it fails with `lazy loading failed for package 'strainhub'`. The one function used, `make_transnet()`, is vendored in [`Tools/strainhub_transnet.R`](Tools/strainhub_transnet.R) (Apache 2.0, with citation and deviations in its header).

---

## Directory Structure

```
VIA/
├── Input_Files/               ← everything you edit (build inputs)
│   ├── genes/                 ← one folder per gene = one entry in the dropdown
│   │   └── COI/
│   │       ├── metadata.csv
│   │       ├── tree.newick
│   │       └── tree.svg
│   ├── MTM_INVASIVE_VECTOR_SPECIES_*.xlsx   ← newest one is used
│   ├── geocode_cache.csv      ← created by the build when it geocodes
│   ├── africa_railways_raw.rds
│   └── shipping_routes_raw.rds
├── Tools/
│   ├── build_app_cache.R      ← the one command
│   └── strainhub_transnet.R
├── Dashboard/                 ← this folder is what gets deployed
│   ├── app.R
│   ├── pubmed_module.R
│   ├── data/                  ← generated by the build; don't edit
│   └── www/                   ← logo, styles; www/genes/ is generated
├── Railway_Extraction/, Shipping_Extraction/, Scripts/, RAXML/   ← upstream steps
```

`Dashboard/` never reads anything outside itself, which is what lets it be deployed on its own.

---

## How to Run

```r
shiny::runApp("Dashboard")
```

Or open `Dashboard/app.R` in RStudio and click **Run App**.

---

## Publishing to Shiny Server / shinyapps.io / Posit Connect

1. Run the build (above).
2. Deploy the **`Dashboard` directory only**:

   ```r
   rsconnect::deployApp("Dashboard", appName = "vector-invasion-atlas")
   ```

   For a self-hosted Shiny Server, copy `Dashboard/` into the app root (e.g. `/srv/shiny-server/vector-invasion-atlas/`) with `data/` and `www/` intact. `sf` needs the GDAL/GEOS/PROJ system libraries on that server.

The server needs outbound internet for PubMed and map tiles.

### Basemap tiles

CARTO now requires an API key for their raster basemap tiles, and serves tiles stamped *"requires an API key"* without one. The app works either way:

| `CARTO_API_KEY` | Base layer | Zoom |
|---|---|---|
| unset (default) | Esri Light Gray Canvas — no key needed | native to 16, upscaled to 19 |
| set | CartoDB Positron | native to 20 |

A key is free at <https://carto.com/basemaps/apikey>. Set it at the top of `app.R`, or pass it at deploy time: `rsconnect::deployApp("Dashboard", envVars = "CARTO_API_KEY")`.

### Performance

Shiny runs each app in a single-threaded R process, so **R time spent on one visitor blocks every other visitor** on that process. The design keeps that near zero:

- All geometry is **serialised to GeoJSON once, at build time**. A visitor's browser receives a ready-made string and parses it itself.
- Railways appear a moment after the map paints; shipping routes are only sent when ticked, and only re-sent when the year filter actually changes which routes are selected.
- Network layouts are precomputed, so the server needs no `igraph`, and the country network looks the same on every visit.

| | Before | Now |
|---|---|---|
| R time per visitor, railways | 2.6 s | **< 0.01 s** |
| R time when shipping is ticked | 0.7 s | **0.04 s** |
| Process startup (libraries + data) | 6.6 s | **2.8 s** |
| Memory per R process | 200 MB | **77 MB** |
| Packages loaded by the app | 26 | **12** |

---

## Using the App

- **Gene** (top of the page) — switches the map, year range, tree, and networks.
- **Pie charts** — one per collection site, split by haplotype. Click one for a breakdown; with Shipping Routes on, that country's routes turn red.
- **Temporal filter** — exact year, up to a year (with playback), a range, or all years.
- **Map layers** — Base / Satellite, plus Land Use (Mapbox), Country Status, Railways and Shipping Routes.
