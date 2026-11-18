# ================================================================================
# build_app_cache.R — rebuild everything the dashboard reads.
#
#   Rscript Tools/build_app_cache.R          (run from the repository root)
#
# Run it whenever you add samples, add a gene, or update the invasive-status
# spreadsheet. It reads only Input_Files/ and writes only Dashboard/data/ and
# Dashboard/www/genes/, so Dashboard/ is ready to deploy as soon as it finishes.
#
# Inputs
#   Input_Files/genes/<GENE>/metadata.csv  one row per sample; its genotype column
#                                          is the one named after the folder
#                                          ("P47Rec Genotype" for P47rec)  (required)
#   Input_Files/genes/<GENE>/tree.newick   RAxML best tree                (optional)
#   Input_Files/genes/<GENE>/tree.svg      iTOL export of that tree       (optional)
#   Input_Files/MTM_INVASIVE_*.xlsx        newest file is used
#   Input_Files/africa_railways_raw.rds    <- Railway_Extraction/prepare_africa_railways.R
#   Input_Files/shipping_routes_raw.rds    <- Shipping_Extraction/prepare_shipping_routes.R
#   Input_Files/Genomics_stats.xlsx       WGS + ddRAD statistics and their descriptions (optional)
#
# Outputs
#   Dashboard/data/genes.rds               per gene: map table, palette, networks
#   Dashboard/data/map_layers.rds          country status, pre-serialised GeoJSON, genomic statistics
#   Dashboard/data/pubmed_snapshot.rds     Literature tab seed (left as-is if offline)
#   Dashboard/www/genes/<GENE>/tree.svg    copied from the gene folder
#   Input_Files/geocode_cache.csv          addresses already looked up — commit it
#
# Every folder under Input_Files/genes/ becomes an entry in the app's gene
# dropdown, named after the folder.
# ================================================================================

# ---- Packages: fail up front with one install line, not halfway through ------
need <- c("sf", "dplyr", "readxl", "rnaturalearth", "rnaturalearthdata", "jsonlite",
          "visNetwork", "igraph", "ape", "tidytree", "castor", "hash", "plyr",
          "tibble", "magrittr", "tidygeocoder", "readr", "shiny", "httr", "xml2",
          "stringr")
absent <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(absent)) {
  stop("Install these packages first:\n  install.packages(c(",
       paste0('"', absent, '"', collapse = ", "), "))", call. = FALSE)
}

suppressPackageStartupMessages({ library(sf); library(dplyr) })
sf::sf_use_s2(FALSE)

if (!dir.exists("Input_Files/genes") || !dir.exists("Dashboard")) {
  stop("Run this from the repository root (the folder containing Input_Files/ and Dashboard/).",
       call. = FALSE)
}

GENES_DIR      <- "Input_Files/genes"
DATA_DIR       <- "Dashboard/data"
WWW_GENES      <- "Dashboard/www/genes"
GEOCODE_CACHE  <- "Input_Files/geocode_cache.csv"
RAW_RAILWAYS   <- "Input_Files/africa_railways_raw.rds"
RAW_SHIPPING   <- "Input_Files/shipping_routes_raw.rds"
GENOMICS_XLSX  <- "Input_Files/Genomics_stats.xlsx"

# Every gene's metadata.csv needs these (either header style works — see
# read_gene()), plus its own genotype/haplotype column, plus FASTA.ID if the
# gene has a tree.newick.
BASE_COLS      <- c("Year", "Country", "Latitude", "Longitude", "City.Town", "County.State")
# Rows whose location matches these are lab colonies, not field sites: kept for
# the networks, left off the map.
INVALID_LOCATIONS <- c("Laboratory")

SIMPLIFY_TOL_M        <- 500   # railway simplification, metres
RAILWAY_MIN_SEGMENT_M <- 500   # drop OSM fragments shorter than this
COORD_DP              <- 4     # ~11 m; shrinks the JSON sent to browsers

fmt <- function(bytes) format(structure(as.numeric(bytes), class = "object_size"), units = "auto")

# ================================================================================
# Geometry helpers
# ================================================================================

drop_degenerate <- function(x) {
  keep <- vapply(sf::st_geometry(x), function(g)
    !sf::st_is_empty(g) && nrow(sf::st_coordinates(g)) > 1, logical(1))
  x[keep, ]
}

# st_simplify() is a no-op on lon/lat (sf won't drop topology on ellipsoidal
# coordinates), so simplify in metres in EPSG:3857 and project back.
simplify_m <- function(x, tol = SIMPLIFY_TOL_M) {
  x <- drop_degenerate(x) %>% sf::st_transform(3857) %>%
    sf::st_simplify(dTolerance = tol, preserveTopology = FALSE)
  sf::st_transform(x[!sf::st_is_empty(x), ], 4326)
}

# One GeoJSON Feature as a string. The app hands these straight to leaflet's
# addGeoJSON(), which JSON.parse()s them in the browser — so R never has to
# serialise this geometry per visitor. That was costing 2.6 s of R time for
# every new session with the railways, blocking everyone else on the process.
feature_json <- function(g, props) {
  geom <- if (inherits(g, "LINESTRING")) {
    list(type = "LineString", coordinates = round(unclass(g)[, 1:2], COORD_DP))
  } else {
    list(type = "MultiLineString",
         coordinates = lapply(unclass(g), function(m) round(m[, 1:2], COORD_DP)))
  }
  if (!length(props)) props <- setNames(list(), character(0))   # {} not []
  as.character(jsonlite::toJSON(list(type = "Feature", properties = props, geometry = geom),
                                auto_unbox = TRUE, digits = NA))
}

feature_collection <- function(features) {
  paste0('{"type":"FeatureCollection","features":[', paste(features, collapse = ","), "]}")
}

# ================================================================================
# Per-gene helpers
# ================================================================================

norm_country <- function(x) tolower(gsub("[^A-Za-z]", "", x))   # "Saudi Arabia" == "SaudiArabia"

# Read one gene's metadata.csv and bring it to a common shape. Stops with a
# plain-language message on anything that would make the dashboard wrong.
read_gene <- function(dir) {
  gene <- basename(dir)
  f <- file.path(dir, "metadata.csv")
  if (!file.exists(f)) stop(sprintf("[%s] missing %s", gene, f), call. = FALSE)

  md <- readr::read_csv(f, show_col_types = FALSE, guess_max = 1e5)
  headers <- names(md)
  # Excel exports name columns "FASTA ID", "City/Town", "Native/Invasive";
  # make.names() turns those into the original COI file's "FASTA.ID",
  # "City.Town", "Native.Invasive", so either style can be dropped in.
  names(md) <- make.names(headers)

  # The gene's own marker is the column named after the gene — "P47Rec Genotype"
  # in folder P47rec — or failing that a plain "Haplotype" column (the original
  # COI file). Other genes' columns, like "COI Haplotype", are left alone.
  hit <- which(startsWith(tolower(names(md)), tolower(make.names(gene))))
  if (length(hit) > 1) {
    stop(sprintf("[%s] more than one column starts with '%s' (%s). Only the genotype column should.",
                 gene, gene, paste(headers[hit], collapse = ", ")), call. = FALSE)
  }
  if (!length(hit)) hit <- which(names(md) == "Haplotype")
  if (!length(hit)) {
    stop(sprintf("[%s] can't find this gene's genotype column. Name it after the folder, e.g. '%s Genotype', or call it 'Haplotype'.",
                 gene, gene), call. = FALSE)
  }
  marker <- names(md)[hit]
  label  <- if (marker == "Haplotype") "Haplotypes" else headers[hit]
  if (marker != "Haplotype") {
    md$Haplotype <- NULL
    names(md)[names(md) == marker] <- "Haplotype"
  }

  has_tree <- file.exists(file.path(dir, "tree.newick"))
  need <- c(BASE_COLS, if (has_tree) "FASTA.ID")
  lost <- setdiff(need, names(md))
  if (length(lost)) {
    stop(sprintf("[%s] metadata.csv is missing column(s): %s\n  Expected at least: %s, and a genotype column",
                 gene, paste(lost, collapse = ", "), paste(need, collapse = ", ")), call. = FALSE)
  }

  # The same specimen entered twice — same Sample ID, genotype, place and year —
  # is counted once. The first P47rec file had 25 of these (Jijiga, 2022: one
  # copy with the COI haplotype filled in, one without), which would otherwise
  # double those pies. The more complete copy is kept. Rows that disagree on
  # anything are left alone and reported by the duplicate check just below.
  if ("Sample.ID" %in% names(md)) {
    ord  <- order(-rowSums(!is.na(md)))   # most complete first
    # Matched on the town, not the coordinates: one copy often has coordinates
    # and the other is blank (SM226 in P47rec), and they are still one specimen.
    key  <- with(md, paste(Sample.ID, Haplotype, Country, City.Town, Year, sep = "\r"))
    twin <- logical(nrow(md))
    twin[ord] <- !is.na(md$Sample.ID[ord]) & !is.na(md$Haplotype[ord]) & duplicated(key[ord])
    if (any(twin)) {
      message(sprintf("  [%s] counted %d sample(s) listed twice only once (same Sample ID, genotype, place and year)",
                      gene, sum(twin)))
      md <- md[!twin, ]
    }
  }

  for (id in intersect(c("Sample.ID", "FASTA.ID"), names(md))) {
    v <- md[[id]][!is.na(md[[id]])]
    dups <- unique(v[duplicated(v)])
    if (length(dups)) warning(sprintf("[%s] duplicate %s(s): %s", gene, id,
                                      paste(head(dups, 10), collapse = ", ")), call. = FALSE)
  }

  # FASTA IDs end in their country ("Kassala8_Sudan"). If that disagrees with the
  # row's Country, the FASTA ID column has slipped out of line with the rest of
  # the sheet — as it had in the first P47rec file. The map and networks don't
  # use FASTA IDs, but a tree is matched on them, so with a tree this is fatal.
  if ("FASTA.ID" %in% names(md)) {
    chk <- !is.na(md$FASTA.ID) & grepl("_", md$FASTA.ID) & !is.na(md$Country)
    bad <- which(chk & norm_country(sub("^.*_", "", md$FASTA.ID)) != norm_country(md$Country))
    if (length(bad)) {
      msg <- sprintf(paste0("[%s] FASTA.ID doesn't match the row's Country on %d of %d rows ",
                            "(e.g. '%s' on a row in %s). The FASTA ID column has probably been ",
                            "sorted or pasted separately from the other columns."),
                     gene, length(bad), sum(chk), md$FASTA.ID[bad[1]], md$Country[bad[1]])
      if (has_tree) stop(msg, " tree.newick is matched on FASTA.ID, so the StrainHub network would be wrong.", call. = FALSE)
      warning(msg, " Not used by the map or networks; it would matter if a tree.newick is added.", call. = FALSE)
    }
  }

  list(gene = gene, dir = dir, md = md, label = label, has_tree = has_tree)
}

# A sample that already has coordinates in any gene's file is at the same place
# in every other gene's file. Using those first means new genes inherit the
# hand-checked coordinates rather than an OpenStreetMap town centroid, and need
# no network calls. Only Sample IDs with a single, unambiguous location are used.
known_locations <- function(all_genes) {
  k <- do.call(rbind, lapply(all_genes, function(g) {
    if (!"Sample.ID" %in% names(g$md)) return(NULL)
    d <- g$md[!is.na(g$md$Sample.ID) & !is.na(g$md$Latitude) & !is.na(g$md$Longitude),
              c("Sample.ID", "Latitude", "Longitude")]
    as.data.frame(d)
  }))
  if (is.null(k)) return(k)
  k <- unique(k)
  k[!k$Sample.ID %in% k$Sample.ID[duplicated(k$Sample.ID)], ]
}

fill_from_known <- function(md, known, gene) {
  if (is.null(known) || !"Sample.ID" %in% names(md)) return(md)
  miss <- which(is.na(md$Latitude) | is.na(md$Longitude))
  hit  <- match(md$Sample.ID[miss], known$Sample.ID)
  ok   <- !is.na(hit)
  if (any(ok)) {
    md$Latitude[miss[ok]]  <- known$Latitude[hit[ok]]
    md$Longitude[miss[ok]] <- known$Longitude[hit[ok]]
    message(sprintf("  [%s] %d sample(s) placed using coordinates already recorded for the same Sample ID",
                    gene, sum(ok)))
  }
  md
}

address_of <- function(d) {
  apply(d[, c("City.Town", "County.State", "Country")], 1, function(r) {
    r <- r[!is.na(r) & nzchar(trimws(r))]
    paste(r, collapse = ", ")
  })
}

# Fill missing coordinates, looking each distinct address up only once, ever.
# The old approach called OpenStreetMap once per row; 540 missing rows here were
# only 25 distinct places. Successful lookups go in geocode_cache.csv so reruns
# make no network calls; failures are not cached, so they retry next time.
geocode_missing <- function(md, gene) {
  miss <- is.na(md$Latitude) | is.na(md$Longitude)
  if (!any(miss)) return(md)

  addr  <- address_of(md[miss, ])
  cache <- if (file.exists(GEOCODE_CACHE)) {
    read.csv(GEOCODE_CACHE, stringsAsFactors = FALSE)
  } else data.frame(address = character(), lat = numeric(), lon = numeric())

  todo <- setdiff(unique(addr[nzchar(addr)]), cache$address)
  if (length(todo)) {
    message(sprintf("  [%s] geocoding %d new place(s) via OpenStreetMap (~1 s each)...",
                    gene, length(todo)))
    got <- tryCatch(tidygeocoder::geo(address = todo, method = "osm", quiet = TRUE),
                    error = function(e) { message("  geocoding failed: ", conditionMessage(e)); NULL })
    if (!is.null(got)) {
      new <- data.frame(address = todo, lat = got$lat, lon = got$long)
      cache <- rbind(cache, new[!is.na(new$lat) & !is.na(new$lon), ])
      write.csv(cache, GEOCODE_CACHE, row.names = FALSE)
    }
  }

  hit <- match(addr, cache$address)
  md$Latitude[miss]  <- cache$lat[hit]
  md$Longitude[miss] <- cache$lon[hit]
  unresolved <- unique(addr[is.na(hit)])
  if (length(unresolved)) {
    warning(sprintf("[%s] %d sample(s) still have no coordinates and are left off the map. ",
                    gene, sum(is.na(hit))),
            "Fix the address or type the coordinates into metadata.csv:\n  ",
            paste(unresolved, collapse = "\n  "), call. = FALSE)
  }
  md
}

haplotype_palette <- function(haps) {
  base <- c("#E41A1C", "#4DAF4A", "#377EB8", "#FF7F00", "#984EA3", "#A65628",
            "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
            "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02",
            "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
  haps <- sort(unique(haps))   # alphabetical (HH, HQ, QQ), so colours don't shift as rows are added
  num  <- suppressWarnings(as.numeric(sub("^\\D*(\\d+).*$", "\\1", haps)))   # Hap12 -> 12
  if (!all(is.na(num))) haps <- haps[order(num, na.last = TRUE)]
  cols <- if (length(haps) <= length(base)) base[seq_along(haps)] else colorRampPalette(base)(length(haps))
  setNames(cols, haps)
}

build_strainhub <- function(md, tree_file) {
  tree <- ape::read.tree(tree_file)
  md <- as.data.frame(md)
  names(md)[names(md) == "FASTA.ID"]     <- "Accession"
  names(md)[names(md) == "Accession.ID"] <- "Unique"
  pruned <- tidytree::drop.tip(tree, setdiff(tree$tip.label, md$Accession))

  graph <- make_transnet(pruned, md, columnSelection = "Country",
                         centralityMetric = 3,   # betweenness
                         treeType = "parsimonious", metricsOutputFile = "")

  nodes <- graph$x$nodes %>% mutate(shape = "dot", color = "grey", font.size = 35)
  # `value` (number of inferred changes) is kept so vis.js scales arrow width by
  # it. The old fixed widths covered only 1-3 changes, so the 25 edges with more
  # had no width and drew thinnest.
  lab <- setNames(nodes$label, nodes$id)
  edges <- graph$x$edges %>%
    mutate(arrows = "to", smooth = TRUE, color = "black", length = 500,
           title  = paste0(lab[as.character(from)], " → ", lab[as.character(to)], ": ",
                           value, " inferred change", ifelse(value == 1, "", "s")))
  list(nodes = nodes, edges = edges)
}

build_country_net <- function(md) {
  clean <- md %>%
    filter(!is.na(Country),   Country   != "", Country   != "NA",
           !is.na(Haplotype), Haplotype != "", Haplotype != "NA")

  set.seed(123)
  sub <- clean %>% group_by(Country, Year) %>% slice_sample(n = 20) %>% ungroup()

  hb <- table(sub$Haplotype, sub$Country); hb[hb > 0] <- 1
  cm <- as.matrix(t(hb) %*% hb); diag(cm) <- 0

  totals <- sub %>% count(Country, name = "n")
  nodes <- data.frame(id = seq_len(nrow(totals)), label = totals$Country,
                      group = totals$Country, value = totals$n,
                      title = paste("Included Samples:", totals$n))
  idx <- which(cm > 0, arr.ind = TRUE)
  idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
  links <- data.frame(from = idx[, 1], to = idx[, 2], value = cm[idx])

  # Bake the igraph layout in here, so the server needs no igraph. The app used
  # visLayout(randomSeed = 1), which seeds vis.js, not igraph — so the layout
  # was different on every R session. randomSeed below fixes it.
  laid <- visNetwork::visNetwork(nodes, links) %>%
    visNetwork::visIgraphLayout(layout = "layout_nicely", randomSeed = 1)
  list(nodes = laid$x$nodes, links = links)
}

# ================================================================================
# 1. GENES
# ================================================================================
message("\n[1/6] Genes")
source("Tools/strainhub_transnet.R")

gene_dirs <- list.dirs(GENES_DIR, recursive = FALSE)
if (!length(gene_dirs)) stop("No gene folders found in ", GENES_DIR, call. = FALSE)

# Read and check every gene before doing any slow work, so a bad file fails fast.
all_genes <- lapply(gene_dirs, read_gene)
known     <- known_locations(all_genes)

unlink(WWW_GENES, recursive = TRUE)   # generated; rebuilt from scratch each run
genes <- list()

for (g in all_genes) {
  gene <- g$gene; dir <- g$dir; md <- g$md

  # ---- Map table ----
  bad <- paste(INVALID_LOCATIONS, collapse = "|")
  is_lab <- grepl(bad, md$City.Town) | grepl(bad, md$County.State)
  map <- md[!is_lab & !is.na(md$Haplotype), ] %>%
    fill_from_known(known, gene) %>%
    geocode_missing(gene)
  map <- map[!is.na(map$Latitude) & !is.na(map$Longitude),
             c("Latitude", "Longitude", "Haplotype", "City.Town", "County.State",
               "Country", "Year")]
  map <- as.data.frame(map)

  ctc <- map[!is.na(map$City.Town) & !is.na(map$Country), c("City.Town", "Country")]
  ctc <- ctc[!duplicated(ctc), ]

  # ---- Networks ----
  tree_file <- file.path(dir, "tree.newick")
  strainhub <- if (file.exists(tree_file)) build_strainhub(md, tree_file) else NULL

  # ---- Tree figure ----
  svg <- file.path(dir, "tree.svg")
  tree_svg <- NULL
  if (file.exists(svg)) {
    dir.create(file.path(WWW_GENES, gene), recursive = TRUE, showWarnings = FALSE)
    file.copy(svg, file.path(WWW_GENES, gene, "tree.svg"), overwrite = TRUE)
    tree_svg <- paste0("genes/", gene, "/tree.svg")
  }

  genes[[gene]] <- list(
    marker_label    = g$label,   # legend title: "Haplotypes", "P47Rec Genotype", ...
    haplo           = map,
    palette         = haplotype_palette(map$Haplotype),
    years           = sort(unique(na.omit(map$Year))),
    city_to_country = setNames(ctc$Country, ctc$City.Town),
    strainhub       = strainhub,
    country_net     = build_country_net(md),
    tree_svg        = tree_svg
  )

  message(sprintf("  %-8s %4d samples on map, coloured by '%s' | StrainHub %s | tree figure %s",
                  gene, nrow(map), g$label,
                  if (is.null(strainhub)) "skipped (no tree.newick)" else
                    sprintf("%d nodes", nrow(strainhub$nodes)),
                  if (is.null(tree_svg)) "none (no tree.svg)" else "yes"))
}

# Order: COI first if present, otherwise alphabetical.
genes <- genes[order(names(genes) != "COI", names(genes))]
saveRDS(genes, file.path(DATA_DIR, "genes.rds"), compress = "xz")

# ================================================================================
# 2. COUNTRY STATUS
# ================================================================================
message("\n[2/6] Country status")
xlsx <- tail(sort(Sys.glob("Input_Files/MTM_INVASIVE_*.xlsx")), 1)
if (!length(xlsx)) stop("No Input_Files/MTM_INVASIVE_*.xlsx found", call. = FALSE)

status <- readxl::read_excel(xlsx, sheet = "Data") %>%
  transmute(name_clean = tolower(COUNTRY_NAME), status_clean = tolower(INVASIVE_STATUS)) %>%
  filter(!is.na(name_clean), !is.na(status_clean)) %>%
  count(name_clean, status_clean) %>%
  slice_max(n, n = 1, by = name_clean, with_ties = FALSE) %>%
  select(-n)

country_status <- rnaturalearth::ne_countries(returnclass = "sf") %>%
  transmute(admin, name_clean = tolower(admin)) %>%
  left_join(status, by = "name_clean") %>%
  select(admin, status_clean)
message(sprintf("  %s -> %d countries with a status", basename(xlsx),
                sum(!is.na(country_status$status_clean))))

# ================================================================================
# 3. RAILWAYS
# ================================================================================
message("\n[3/6] Railways")
if (!file.exists(RAW_RAILWAYS)) stop(RAW_RAILWAYS, " not found. Run Railway_Extraction/prepare_africa_railways.R first.", call. = FALSE)

rail <- readRDS(RAW_RAILWAYS) %>% select(fclass)

# OSM splits a line at every attribute change: the median segment is 355 m and
# a quarter are under 100 m (sidings, yards) — invisible at this map's zoom, but
# each costs at least two vertices. 500 m keeps ~95% of the network by length.
seg_m <- as.numeric(sf::st_length(sf::st_transform(rail, 3857)))
rail  <- rail[seg_m >= RAILWAY_MIN_SEGMENT_M, ] %>% simplify_m() %>%
  group_by(fclass) %>%
  # st_combine, not summarise()'s default st_union: union nodes every crossing
  # and added ~50,000 vertices for nothing.
  summarise(geometry = sf::st_combine(geometry), .groups = "drop")

RAIL_COLOURS <- c(rail = "#4a4a4a", light_rail = "#8a6db5",
                  narrow_gauge = "#c0622a", preserved = "#2e7d52")
rail_features <- vapply(seq_len(nrow(rail)), function(i) {
  ft  <- rail$fclass[i]
  col <- if (ft %in% names(RAIL_COLOURS)) RAIL_COLOURS[[ft]] else "#888888"
  feature_json(sf::st_geometry(rail)[[i]],
               list(style = list(color = col), popup = paste("Railway type:", ft)))
}, character(1))
railways_geojson <- feature_collection(rail_features)
message(sprintf("  %d railway types, %s of GeoJSON", nrow(rail), fmt(nchar(railways_geojson))))

# ================================================================================
# 4. SHIPPING ROUTES
# ================================================================================
message("\n[4/6] Shipping routes")
if (!file.exists(RAW_SHIPPING)) stop(RAW_SHIPPING, " not found. Run Shipping_Extraction/prepare_shipping_routes.R first.", call. = FALSE)

ship_raw <- readRDS(RAW_SHIPPING)
flat <- do.call(rbind, lapply(names(ship_raw), function(y) {
  s <- ship_raw[[y]]; if (is.null(s) || !nrow(s)) return(NULL); s$year <- as.integer(y); s
}))
flat$route_id <- paste(flat$from, flat$to, sep = "\r")

# The raw data repeats each from-to pair once per year (34,735 rows, 1,774
# distinct). Store each geometry once plus a route/year index. Not simplified:
# searoute paths hug coastlines and simplifying can walk a lane across land.
routes <- drop_degenerate(flat[!duplicated(flat$route_id), c("route_id", "from", "to")])
ship_features <- setNames(
  vapply(sf::st_geometry(routes), feature_json, character(1), props = list()),
  routes$route_id)

shipping <- list(
  features    = ship_features,
  routes      = sf::st_drop_geometry(routes),
  route_years = flat %>% sf::st_drop_geometry() %>% distinct(route_id, year) %>%
                  filter(route_id %in% routes$route_id),
  years       = sort(unique(flat$year))
)
message(sprintf("  %d unique routes (from %d rows), %s of GeoJSON", nrow(routes),
                nrow(flat), fmt(sum(nchar(ship_features)))))

# ================================================================================
# 5. GENOMIC STATISTICS (optional) — WGS and ddRAD boxes, map pop-ups, About
# ================================================================================
# Sheets of GENOMICS_XLSX:
#   "Data Descriptors"  Statistic | Scale | Tool Tip | Info Panel, then a
#                       "How to read these genomic statistics" row
#   "WGS"               one row per country
#   "ddRAD"             per-site table; below a blank row, a one-line note and a
#                       per-region table whose header row starts with "Region"
# Columns are found by name pattern, so fixing a header's spelling doesn't break it.
message("\n[5/6] Genomic statistics")
genomics <- NULL
if (file.exists(GENOMICS_XLSX)) {
  lost <- setdiff(c("Data Descriptors", "WGS", "ddRAD"), readxl::excel_sheets(GENOMICS_XLSX))
  if (length(lost)) stop(GENOMICS_XLSX, " is missing sheet(s): ", paste(lost, collapse = ", "), call. = FALSE)

  # The one column of `d` whose header matches `pattern`, as numbers if `num`.
  pick <- function(d, pattern, sheet, num = FALSE) {
    i <- grep(pattern, names(d), ignore.case = TRUE, perl = TRUE)
    if (length(i) != 1) stop(GENOMICS_XLSX, ", sheet ", sheet, ": expected one column matching '",
                             pattern, "', found ", length(i), call. = FALSE)
    x <- d[[i]]
    if (num) suppressWarnings(as.numeric(x)) else trimws(as.character(x))
  }
  # "Samake et al. 2023" (+ " (dissertation)"), one per ";"-separated citation.
  short_cite <- function(src) vapply(strsplit(src, ";\\s*(?=[A-Z])", perl = TRUE), function(p) {
    paste(paste0(sub("^(\\w+).*$", "\\1", p),
                 ifelse(grepl("&|et al", p), " et al.", ""), " ",
                 sub("^.*?\\b((19|20)\\d{2})\\b.*$|^.*$", "\\1", p, perl = TRUE),
                 ifelse(grepl("Order No|thesis|dissertation", p, ignore.case = TRUE), " (dissertation)", "")),
          collapse = "; ")
  }, character(1))
  first_doi <- function(src) sub("\\.$", "", ifelse(grepl("https?://doi\\.org/", src),
                                  sub("^.*?(https?://doi\\.org/\\S+).*$", "\\1", src), NA_character_))
  # Countries spelled as in the gene metadata, for highlighting the clicked pie's country.
  known_countries <- unique(unlist(lapply(all_genes, function(g) g$md$Country)))
  check_countries <- function(x, sheet) {
    bad <- unique(x[!norm_country(x) %in% norm_country(known_countries)])
    if (length(bad)) warning(GENOMICS_XLSX, ", sheet ", sheet, ": ", paste0("'", bad, "'", collapse = ", "),
                             " not found in any gene's metadata.csv; check the spelling (its row ",
                             "can't be highlighted on the map).", call. = FALSE)
  }

  # ---- Descriptors ----
  dd <- as.data.frame(readxl::read_excel(GENOMICS_XLSX, sheet = "Data Descriptors", .name_repair = "unique_quiet"))
  stat <- pick(dd, "^Statistic", "Data Descriptors")
  tip  <- pick(dd, "Tool ?Tip", "Data Descriptors")
  info <- pick(dd, "Info ?Panel", "Data Descriptors")
  key  <- ifelse(grepl("Tajima", stat), "tajima_d", ifelse(grepl("π|Nucleotide", stat), "pi",
                 ifelse(grepl("FST", stat, ignore.case = TRUE), "fst", NA)))
  descriptors <- data.frame(key = key, label = stat, tip = tip, info = info)[!is.na(key) & !is.na(tip), ]
  how <- which(grepl("^How to read", stat))
  reading <- if (length(how)) unname(na.omit(unlist(dd[how[1], -1]))[1]) else NA_character_
  if (!setequal(descriptors$key, c("tajima_d", "pi", "fst")))
    stop(GENOMICS_XLSX, ", sheet Data Descriptors: expected rows for Tajima's D, π and FST", call. = FALSE)

  # ---- WGS ----
  w <- as.data.frame(readxl::read_excel(GENOMICS_XLSX, sheet = "WGS", .name_repair = "unique_quiet"))
  wgs <- data.frame(country  = pick(w, "^(Population|Country)", "WGS"),
                    sites    = pick(w, "Sites", "WGS", TRUE),
                    n        = pick(w, "Sample", "WGS", TRUE),
                    tajima_d = pick(w, "Tajima", "WGS", TRUE),
                    pi       = pick(w, "Diversity|π", "WGS", TRUE),
                    fst      = pick(w, "FST", "WGS", TRUE),
                    source   = pick(w, "Source", "WGS"))
  wgs <- wgs[!is.na(wgs$country) & wgs$country != "", ]
  check_countries(wgs$country, "WGS")

  # ---- ddRAD: site table, note, region table, stacked on one sheet ----
  g <- as.data.frame(readxl::read_excel(GENOMICS_XLSX, sheet = "ddRAD", col_names = FALSE,
                                        .name_repair = "unique_quiet"))
  g[] <- lapply(g, function(x) trimws(as.character(x)))
  blank <- rowSums(!is.na(g)) == 0
  block <- function(from) {           # header row `from` to the next blank row
    to <- c(which(blank & seq_along(blank) > from), nrow(g) + 1)[1] - 1
    d  <- g[seq(from + 1, length.out = max(0, to - from)), , drop = FALSE]
    names(d) <- unlist(g[from, ])
    d[, !is.na(names(d)), drop = FALSE]
  }
  s <- block(1)
  ddrad <- data.frame(site     = pick(s, "^Site$", "ddRAD"),
                      country  = pick(s, "^Country$", "ddRAD"),
                      lat      = pick(s, "^Lat", "ddRAD", TRUE),
                      lon      = pick(s, "^Lon", "ddRAD", TRUE),
                      n        = pick(s, "Sample", "ddRAD", TRUE),
                      pi       = pick(s, "^(π|pi)$|Diversity", "ddRAD", TRUE),
                      S        = pick(s, "^S$", "ddRAD", TRUE),
                      tajima_d = pick(s, "Tajima", "ddRAD", TRUE),
                      source   = pick(s, "Source", "ddRAD"))
  ddrad <- ddrad[!is.na(ddrad$site) & !is.na(ddrad$lat) & !is.na(ddrad$lon), ]
  ddrad$doi   <- first_doi(ddrad$source)
  ddrad$short <- short_cite(ddrad$source)
  check_countries(ddrad$country, "ddRAD")

  hdr <- which(g[[1]] == "Region")[1]
  row  <- seq_len(nrow(g))
  note <- na.omit(g[[1]][row > nrow(s) + 1 & (is.na(hdr) | row < hdr)])
  regions <- NULL
  if (!is.na(hdr)) {
    r <- block(hdr)
    regions <- data.frame(region   = pick(r, "^Region$", "ddRAD"),
                          sites    = pick(r, "Sites", "ddRAD", TRUE),
                          n        = pick(r, "Sample", "ddRAD", TRUE),
                          pi       = pick(r, "π|Diversity", "ddRAD", TRUE),
                          tajima_d = pick(r, "Tajima", "ddRAD", TRUE),
                          fst      = pick(r, "FST", "ddRAD", TRUE),
                          source   = pick(r, "Source", "ddRAD"))
    regions <- regions[!is.na(regions$region), ]
    # "NE Ethiopia" -> Ethiopia: the site-table country named in the region.
    cs <- unique(ddrad$country)
    regions$country <- vapply(regions$region, function(x) c(cs[vapply(cs, grepl, logical(1), x = x, fixed = TRUE)], NA)[1], "")
    if (anyNA(regions$country)) warning(GENOMICS_XLSX, ", sheet ddRAD: no country found in region(s) ",
                                        paste(regions$region[is.na(regions$country)], collapse = ", "), call. = FALSE)
    regions$short <- short_cite(regions$source)
  }

  genomics <- list(descriptors = descriptors, reading = reading, wgs = wgs,
                   ddrad_sites = ddrad, ddrad_regions = regions, ddrad_note = paste(note, collapse = " "))
  message(sprintf("  WGS: %d countries | ddRAD: %d sites, %d regions | %d descriptors",
                  nrow(wgs), nrow(ddrad), NROW(regions), nrow(descriptors)))
} else {
  message("  none (", GENOMICS_XLSX, " not found)")
}

# ---- Self-check: everything the browser will JSON.parse() actually parses ----
stopifnot(
  !is.null(jsonlite::fromJSON(railways_geojson, simplifyVector = FALSE)$features),
  length(jsonlite::fromJSON(feature_collection(ship_features),
                            simplifyVector = FALSE)$features) == nrow(routes),
  all(vapply(genes, function(g) all(g$haplo$Haplotype %in% names(g$palette)), logical(1))),
  all(vapply(genes, function(g) !anyNA(g$haplo[, c("Latitude", "Longitude")]), logical(1)))
)

saveRDS(list(country_status   = country_status,
             railways_geojson = railways_geojson,
             railway_types    = sort(unique(rail$fclass)),
             railway_colours  = RAIL_COLOURS,
             shipping         = shipping,
             genomics         = genomics),
        file.path(DATA_DIR, "map_layers.rds"), compress = "xz")

# ================================================================================
# 6. PUBMED SNAPSHOT
# ================================================================================
message("\n[6/6] PubMed snapshot")
old <- setwd("Dashboard"); source("pubmed_module.R"); setwd(old)
papers <- tryCatch(.get_pubmed_papers(force = TRUE), error = function(e) NULL)
if (!is.null(papers) && nrow(papers) > 0) {
  saveRDS(list(data = papers, fetched = Sys.time()),
          file.path(DATA_DIR, "pubmed_snapshot.rds"), compress = "xz")
  message(sprintf("  %d papers", nrow(papers)))
} else {
  message("  PubMed unreachable — keeping the previous snapshot")
}

# ================================================================================
message("\nDone. Dashboard/ is ready to deploy:")
for (f in c("genes.rds", "map_layers.rds", "pubmed_snapshot.rds")) {
  p <- file.path(DATA_DIR, f)
  if (file.exists(p)) message(sprintf("  data/%-20s %s", f, fmt(file.info(p)$size)))
}
