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

### Update the genomic statistics (WGS and ddRAD)

Edit `Input_Files/Genomics_stats.xlsx` and run the build. It needs three sheets:

| Sheet | What it drives |
|---|---|
| `Data Descriptors` | Hover text (`Tool Tip`) and the "How to read these statistics" panel (`Info Panel`) for Tajima's D, π and FST, plus the "How to read these genomic statistics" row |
| `WGS` | The map sidebar's **WGS data** table, one row per country |
| `ddRAD` | Per-site table (map pop-ups); below a blank row, a one-line note and a per-region table whose header starts with `Region` (the **ddRAD data** table) |

Columns are found by name (`Tajima`, `FST`, `Sample`…), so fixing a header's spelling is safe. Country names must match the `Country` column of the gene metadata (the build warns if one doesn't), or that row can't be highlighted when a pie chart is clicked. The WGS values come from `Scripts/Genomic_population_statistics.sh`.

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
- **WGS population statistics** — `Scripts/Genomic_population_statistics.sh` (vcftools on a PBS cluster).
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
│   ├── Genomics_stats.xlsx    ← WGS + ddRAD statistics (optional)
│   ├── geocode_cache.csv      ← created by the build when it geocodes
│   ├── africa_railways_raw.rds
│   └── shipping_routes_raw.rds
├── Tools/
│   ├── build_app_cache.R      ← the one command
│   └── strainhub_transnet.R
├── Dashboard/                 ← this folder is what gets deployed
│   ├── app.R
│   ├── pubmed_module.R
│   ├── manifest.json          ← R version + packages for Connect Cloud
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

## Publishing (Posit Connect Cloud)

VIA is hosted on [Posit Connect Cloud](https://connect.posit.cloud), which publishes straight from this GitHub repository. Connect Cloud replaces shinyapps.io. It reads `Dashboard/manifest.json` to know which R version and packages to install, and serves everything in `Dashboard/` (minus what `Dashboard/.rscignore` lists).

### One-time setup (already done unless you are moving accounts)

1. Sign in to <https://connect.posit.cloud> with the **lab's shared account**, not a personal one, so the next student can take over by getting that login. Signing in with the GitHub account that owns this repository is simplest.
2. When Connect Cloud asks, let its GitHub app see this repository.
3. Click **Publish** → **Shiny**, then choose:
   - repository: this one (`sayalg/VIA`), branch `main`
   - primary file: **`Dashboard/app.R`**
4. Turn on **Automatically publish on push**, so every push to `main` republishes the app.
5. Optional: under the content's **Variables**, add `CARTO_API_KEY` (see *Basemap tiles* below).
6. Make sure the content's access is set to **public**, then use its URL everywhere the dashboard is linked (manuscript, About page, this README).

### Updating the live dashboard

| You changed… | Do this, then commit and push |
|---|---|
| Data in `Input_Files/` | `Rscript Tools/build_app_cache.R`, and commit `Dashboard/data/` and `Dashboard/www/genes/` |
| App code only (`app.R`, `pubmed_module.R`, `www/styles.css`) | nothing extra |
| Added a `library()` to the app, or updated R packages | `rsconnect::writeManifest("Dashboard", appPrimaryDoc = "app.R")` and commit `Dashboard/manifest.json` |

Nothing goes live until it's pushed to `main`. If a publish fails, open the content in Connect Cloud and read its **Logs**: a missing package almost always means `manifest.json` is out of date.

On the free plan the app sleeps when nobody is using it, so the first visitor after a quiet spell waits a little while it starts.

### Self-hosting instead

Copy `Dashboard/` to a Shiny Server app folder (e.g. `/srv/shiny-server/vector-invasion-atlas/`) with `data/` and `www/` intact, and install the packages listed under *Required R Packages*. `sf` needs the GDAL/GEOS/PROJ system libraries on that server.

Wherever it runs, the server needs outbound internet for PubMed and map tiles.

### Basemap tiles

CARTO now requires an API key for their raster basemap tiles, and serves tiles stamped *"requires an API key"* without one. The app works either way:

| `CARTO_API_KEY` | Base layer | Zoom |
|---|---|---|
| unset (default) | Esri Light Gray Canvas — no key needed | native to 16, upscaled to 19 |
| set | CartoDB Positron | native to 20 |

A key is free at <https://carto.com/basemaps/apikey>. On Connect Cloud, add it as a `CARTO_API_KEY` variable in the content's settings rather than writing it into `app.R`, so it stays out of the repository.

---

## Using the App

- **Gene** (top of the page) — switches the map, year range, tree, and networks.
- **Pie charts** — one per collection site, split by haplotype. Click one for a breakdown; its country is highlighted in the sidebar's WGS and ddRAD tables, and with Shipping Routes on, that country's routes turn red.
- **Phylogenetic tree** — scroll or use + / − to zoom, drag to pan.
- **The VIA logo/name** — returns to the map from any tab.
- **Temporal filter** — exact year, up to a year (with playback), a range, or all years.
- **Map layers** — Base / Satellite, plus Land Use (Mapbox), Country Status, Railways and Shipping Routes.
