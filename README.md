# VIA: Vector Invasion Atlas

An interactive R Shiny dashboard showing *Anopheles stephensi* haplotype distribution on a map (pie chart per collection site), phylogenetic trees, phylogeographic dispersal and haplotype-sharing networks, and a live hub of *An. stephensi* literature.

Web portal to VIA can be found here: https://carterlab-via.share.connect.posit.cloud/

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

1. Export the gene's sample sheet from Excel as CSV, saved as `metadata.csv`. Use the same columns as the other genes' sheets, plus the gene's own result column, named after the gene: **`<Gene> Genotype`** (e.g. `P47Rec Genotype`) or just **`<Gene>`** (COI's column is `COI`).
2. Put it in a folder named after the gene: `Input_Files/genes/P47rec/metadata.csv`. **The folder name must match the start of that column** (capitals don't matter), because that's how the build knows which column to colour the map by. Any other extra column is carried along and doesn't interfere.
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
| `Data Descriptors` | Hover text (`Tool Tip`) for Tajima's D, π and Fst. The longer descriptions are written in `Dashboard/app.R` (`GENO_TEXT`), not in this sheet |
| `WGS` | The **WGS** map layer, one row per country, drawn as one circle per country |
| `ddRAD` | The **ddRAD** map layer: a per-site table (with coordinates); below a blank row, a one-line note and a per-region table whose header starts with `Region` (shown in each site's pop-up) |

Columns are found by name (`Tajima`, `FST`, `Sample`…), so fixing a header's spelling is safe. WGS country names must match the `Country` column of the gene metadata and a country name on the map (e.g. `Iran`, `Saudi Arabia`), or that country gets no circle. Each ddRAD site is assigned to a region by the region table's site counts, taken in sheet order; the build warns if a region's sample sizes don't add up. The WGS values come from `Scripts/Genomic_population_statistics.sh`.

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
- **Haplotype diversity and sampling sensitivity** — `Scripts/Hapdiversitycalculation.R` (run from the repository root; reads `Input_Files/genes/COI/metadata.csv`, writes to `via_sensitivity_output/`: CSVs of Hd per country, edge recovery across subsamples and rarefied Hd, plus the country haplotype-sharing network as an interactive HTML page with edge stability — dashed red edges appeared in fewer than 95% of 1,000 subsamples). Needs the `visNetwork` and `igraph` packages, and pandoc to save the HTML (RStudio includes it).
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
├── Dockerfile, .dockerignore  ← container image of Dashboard/
├── Dashboard/                 ← this folder is what gets deployed
│   ├── app.R
│   ├── pubmed_module.R
│   ├── manifest.json          ← R version + packages for Connect Cloud
│   ├── data/                  ← generated by the build; don't edit
│   └── www/                   ← logo, styles; www/genes/ is generated
├── Supplementary/             ← supplementary figures and tables for the manuscript
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

### Run it in a container (Docker)

The `Dockerfile` at the repository root packages the dashboard with R 4.5.2 and the packages it needs, so it runs the same on any machine with Docker:

```bash
docker build -t via .
docker run --rm -p 3838:3838 via
```

Then open <http://localhost:3838>. To use the CARTO basemap, add `-e CARTO_API_KEY=your_key` to `docker run`. The image contains only `Dashboard/` and the data already built into it, so after changing `Input_Files/` run `Rscript Tools/build_app_cache.R` and rebuild the image. The first build takes several minutes. The container needs outbound internet for PubMed and map tiles.

### Self-hosting instead

Copy `Dashboard/` to a Shiny Server app folder (e.g. `/srv/shiny-server/vector-invasion-atlas/`) with `data/` and `www/` intact, and install the packages listed under *Required R Packages*. `sf` needs the GDAL/GEOS/PROJ system libraries on that server.

Wherever it runs, the server needs outbound internet for PubMed and map tiles.

## Using the App

- **Genetic Diversity** menu (top of the page) — chooses which gene's data (COI by default) drives the map, year range, tree, and networks.
- **Dispersal Network** and **Haplotype Network** — the phylogeographic dispersal network (StrainHub) and the country haplotype-sharing network, on separate tabs. In the haplotype network, dashed red edges were found in fewer than 95% of 1,000 repeated subsamples.
- **Pie charts** — one per collection site, split by haplotype. Click one for a breakdown; with Shipping Routes on, that country's routes turn red.
- **Show on map** (above the map) — Haplotype (the pie charts), WGS or ddRAD. The WGS and ddRAD layers draw a circle per location, coloured yellow→red by nucleotide diversity (π); click one for its π, Tajima's D, Fst and sample counts. The temporal filter applies to haplotypes only.
- **Phylogenetic tree** — scroll or use + / − to zoom, drag to pan.
- **The VIA logo/name** — returns to the map from any tab.
- **Temporal filter** — exact year, up to a year (with playback), a range, or all years.
- **Map layers** — Base / Satellite, plus Land Use (Mapbox), Country Status, Railways and Shipping Routes.
