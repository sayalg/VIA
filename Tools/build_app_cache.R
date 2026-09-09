# ================================================================================
# build_app_cache.R
#
# Run ONCE offline (and again whenever the underlying data changes). Turns the
# raw pipeline outputs into the small, pre-chewed caches that Dashboard/app.R
# loads at startup, so the deployed app does no heavy computation and — far more
# importantly — ships two orders of magnitude less geometry to each browser.
#
#   Usage, from the repository root:
#       Rscript Tools/build_app_cache.R
#
# Inputs   (raw, reproducible pipeline artefacts)
#   Input_Files/africa_railways_raw.rds     <- Railway_Extraction/prepare_africa_railways.R
#   Input_Files/shipping_routes_raw.rds     <- Shipping_Extraction/prepare_shipping_routes.R
#   Dashboard/data/FinalCOI_metadata_w_coor.csv
#   Dashboard/data/COI_subsampled20_per_country_year.aligned.raxml.bestTree
#
# Outputs  (app caches — these are what get deployed)
#   Dashboard/data/africa_railways.rds      merged by type, simplified, 4 dp
#   Dashboard/data/shipping_routes.rds      deduplicated routes + year index
#   Dashboard/data/networks.rds             both visNetwork node/edge tables
#   Dashboard/data/pubmed_snapshot.rds      seed results for the Literature tab
#
# WHY THIS EXISTS — measured on the raw data:
#   * Shipping "All Years" (the default view) was 34,735 features / 57.2 MB of
#     GeoJSON per browser. 32,959 of those were the *same* from-to pair repeated
#     once per year, drawn on top of itself. Deduplicating is visually lossless
#     and gets it to 1,776 features / 2.5 MB.
#   * Railways were 32,564 separate features / 16 MB. That is 32k SVG paths for
#     Leaflet to lay out, which is what made panning crawl. Merging to one
#     MULTILINESTRING per type gives 3 paths / 2.6 MB.
#   * st_simplify() is a no-op on lon/lat data (sf refuses to drop topology on
#     ellipsoidal coordinates), so geometry is projected to EPSG:3857, simplified
#     with a tolerance in metres, and projected back.
# ================================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
})

sf::sf_use_s2(FALSE)

RAW_RAILWAYS  <- "Input_Files/africa_railways_raw.rds"
RAW_SHIPPING  <- "Input_Files/shipping_routes_raw.rds"
DATA_DIR      <- "Dashboard/data"
METADATA_CSV  <- file.path(DATA_DIR, "FinalCOI_metadata_w_coor.csv")
TREE_FILE     <- file.path(DATA_DIR,
                           "COI_subsampled20_per_country_year.aligned.raxml.bestTree")

# Simplification tolerance in metres. 500 m is invisible at continental zoom,
# where this layer is actually looked at.
SIMPLIFY_TOL_M <- 500

# Railway segments shorter than this are dropped from the cache. See the railway
# section below for why this is the single biggest lever on how long the layer
# takes to appear.
RAILWAY_MIN_SEGMENT_M <- 500

# Coordinate decimal places kept in the cache. 4 dp is ~11 m at the equator.
COORD_DP <- 4

if (!dir.exists(DATA_DIR)) {
  stop("Run this from the repository root — ", DATA_DIR, " not found.")
}

fmt <- function(bytes) {
  format(structure(as.numeric(bytes), class = "object_size"), units = "auto")
}

# Rounding coordinates shrinks the JSON that leaflet serialises to the browser.
# A LINESTRING unclasses to a bare matrix; a MULTILINESTRING to a list of them,
# so the two cases need handling separately.
round_coords <- function(x, dp = COORD_DP) {
  rounded <- lapply(sf::st_geometry(x), function(g) {
    cls <- class(g)
    raw <- unclass(g)
    if (is.matrix(raw)) {
      structure(round(raw, dp), class = cls)
    } else if (is.list(raw)) {
      structure(lapply(raw, function(m) round(m, dp)), class = cls)
    } else {
      g
    }
  })
  sf::st_geometry(x) <- sf::st_sfc(rounded, crs = sf::st_crs(x))
  x
}

# A handful of routes come back from searoute as degenerate single-point
# "lines". GEOS throws on those, so drop them before any geometry op.
drop_degenerate <- function(x) {
  keep <- vapply(sf::st_geometry(x), function(g) {
    !sf::st_is_empty(g) && nrow(sf::st_coordinates(g)) > 1
  }, logical(1))
  x[keep, ]
}

simplify_m <- function(x, tol = SIMPLIFY_TOL_M) {
  x %>%
    drop_degenerate() %>%
    sf::st_transform(3857) %>%
    sf::st_simplify(dTolerance = tol, preserveTopology = FALSE) %>%
    (\(g) g[!sf::st_is_empty(g), ])() %>%
    sf::st_transform(4326)
}

# ================================================================================
# 1. RAILWAYS
# ================================================================================
message("\n[1/4] Railways")

if (!file.exists(RAW_RAILWAYS)) {
  stop(RAW_RAILWAYS, " not found. Run Railway_Extraction/prepare_africa_railways.R first.")
}

rail_raw <- readRDS(RAW_RAILWAYS)
message(sprintf("  raw: %d features, %s in memory",
                nrow(rail_raw), fmt(object.size(rail_raw))))

# app.R colours by fclass and nothing else; every other column is dead weight
# that would be serialised to the browser.
rail_slim <- rail_raw %>% select(fclass)

# Drop very short segments. OpenStreetMap splits a railway at every change of
# attribute, so the raw extract is dominated by fragments: the median segment is
# 355 m and a quarter of them are under 100 m — sidings, yard tracks and
# bookkeeping splits, none of which are visible at the zooms this map is used at.
# They are, however, most of the cost: each one keeps at least two vertices no
# matter how hard the geometry is simplified.
#
# At 500 m the layer keeps ~95% of the network by length for ~40% of the
# vertices. Raise this to trade more of the network away for more speed, or set
# it to 0 to keep everything.
seg_len_m <- as.numeric(sf::st_length(sf::st_transform(rail_slim, 3857)))
kept      <- rail_slim[seg_len_m >= RAILWAY_MIN_SEGMENT_M, ]
message(sprintf("  segments >= %d m: %d of %d (%.1f%% of network length)",
                RAILWAY_MIN_SEGMENT_M, nrow(kept), nrow(rail_slim),
                100 * sum(seg_len_m[seg_len_m >= RAILWAY_MIN_SEGMENT_M]) / sum(seg_len_m)))

africa_railways <- kept %>%
  simplify_m() %>%
  group_by(fclass) %>%
  # st_combine, NOT the st_union that summarise() defaults to. st_union nodes the
  # lines at every intersection, which added ~50,000 vertices here for no visible
  # benefit; st_combine just collects the geometries.
  summarise(geometry = sf::st_combine(geometry), .groups = "drop") %>%
  round_coords()

message(sprintf("  cached: %d features (one per railway type), %d vertices, %s in memory",
                nrow(africa_railways),
                sum(vapply(sf::st_geometry(africa_railways),
                           function(g) nrow(sf::st_coordinates(g)), numeric(1))),
                fmt(object.size(africa_railways))))
saveRDS(africa_railways, file.path(DATA_DIR, "africa_railways.rds"), compress = "xz")

# ================================================================================
# 2. SHIPPING ROUTES
# ================================================================================
message("\n[2/4] Shipping routes")

if (!file.exists(RAW_SHIPPING)) {
  stop(RAW_SHIPPING, " not found. Run Shipping_Extraction/prepare_shipping_routes.R first.")
}

ship_raw <- readRDS(RAW_SHIPPING)
years    <- as.integer(names(ship_raw))
message(sprintf("  raw: %d years, %d total features",
                length(ship_raw), sum(vapply(ship_raw, nrow, integer(1)))))

# Flatten to one table, tagging each row with the year it came from.
flat <- do.call(rbind, lapply(names(ship_raw), function(y) {
  s <- ship_raw[[y]]
  if (is.null(s) || nrow(s) == 0) return(NULL)
  s$year <- as.integer(y)
  s
}))

flat$route_id <- paste(flat$from, flat$to, sep = "\r")

# The year index: which routes were active in which year. Cheap to filter on.
route_years <- flat %>%
  sf::st_drop_geometry() %>%
  distinct(route_id, year)

# One geometry per unique from-to pair. Drawing the same line once per year was
# purely redundant — the extra copies rendered on top of each other.
#
# Deliberately NOT simplified: these paths are routed around coastlines by
# searoute, and shaving vertices off a sea lane can walk it across land. The
# vertex count is already modest (~66 per route), so deduplication and
# coordinate rounding do essentially all the work here anyway.
routes <- flat[!duplicated(flat$route_id), c("route_id", "from", "to")] %>%
  drop_degenerate() %>%
  round_coords()

# Keep the year index consistent with whatever geometry survived.
route_years <- route_years %>% filter(route_id %in% routes$route_id)

message(sprintf("  cached: %d unique routes (was %d), %s in memory",
                nrow(routes), nrow(flat), fmt(object.size(routes))))

saveRDS(list(routes = routes, route_years = route_years, years = sort(unique(years))),
        file.path(DATA_DIR, "shipping_routes.rds"), compress = "xz")

# ================================================================================
# 3. NETWORKS
# ================================================================================
message("\n[3/4] Networks")

source("Dashboard/strainhub_transnet.R")

metadata_full <- readr::read_csv(METADATA_CSV, col_names = TRUE, show_col_types = FALSE)

# ---- Panel 1: StrainHub country transmission network ----
build_haplo <- function() {
  treedata <- ape::read.tree(TREE_FILE)
  metadata <- metadata_full

  names(metadata)[names(metadata) == "FASTA.ID"]     <- "Accession"
  names(metadata)[names(metadata) == "Accession.ID"] <- "Unique"

  taxa_without_data <- setdiff(treedata$tip.label, metadata$Accession)
  pruned_tree <- tidytree::drop.tip(treedata, taxa_without_data)

  metadata$haplo_concat <- paste(metadata$Country, metadata$Haplotype, sep = "_")

  graph <- make_transnet(pruned_tree,
                         metadata,
                         columnSelection   = "Country",
                         centralityMetric  = 3,   # Betweenness Centrality
                         treeType          = "parsimonious",
                         metricsOutputFile = "")

  nodes <- graph$x$nodes %>%
    mutate(shape = "dot", color = "grey", font.size = 35)

  edges <- graph$x$edges %>%
    left_join(nodes %>% select(id, group), by = c("from" = "id")) %>%
    rename(group_from = group) %>%
    left_join(nodes %>% select(id, group), by = c("to" = "id")) %>%
    rename(group_to = group) %>%
    mutate(arrows = "to",
           smooth = TRUE,
           color  = "black",
           width  = case_when(value == 1 ~ 1, value == 2 ~ 5, value == 3 ~ 20),
           value  = NULL,
           length = 500)

  list(nodes = nodes, edges = edges)
}

# ---- Panel 2: Country connectivity (shared haplotypes) ----
build_country <- function() {
  metadata_clean_cn <- metadata_full %>%
    filter(!is.na(Country)   & Country   != "" & Country   != "NA") %>%
    filter(!is.na(Haplotype) & Haplotype != "" & Haplotype != "NA")

  set.seed(123)
  metadata_subsampled_cn <- metadata_clean_cn %>%
    group_by(Country, Year) %>%
    slice_sample(n = 20) %>%
    ungroup()

  h_matrix_raw    <- table(metadata_subsampled_cn$Haplotype, metadata_subsampled_cn$Country)
  h_matrix_binary <- h_matrix_raw
  h_matrix_binary[h_matrix_binary > 0] <- 1

  country_matrix <- as.matrix(t(h_matrix_binary) %*% h_matrix_binary)
  diag(country_matrix) <- 0

  country_totals_cn <- metadata_subsampled_cn %>%
    group_by(Country) %>%
    summarize(n = n(), .groups = "drop")

  nodes_country <- data.frame(
    id    = seq_len(nrow(country_totals_cn)),
    label = country_totals_cn$Country,
    group = country_totals_cn$Country,
    value = country_totals_cn$n,
    title = paste("Included Samples:", country_totals_cn$n)
  )

  edges_idx <- which(country_matrix > 0, arr.ind = TRUE)
  edges_idx <- edges_idx[edges_idx[, 1] < edges_idx[, 2], , drop = FALSE]
  links_country <- data.frame(
    from  = edges_idx[, 1],
    to    = edges_idx[, 2],
    value = country_matrix[edges_idx]
  )

  list(nodes = nodes_country, links = links_country)
}

haplo_net   <- build_haplo()
country_net <- build_country()

message(sprintf("  strainhub network : %d nodes, %d edges",
                nrow(haplo_net$nodes), nrow(haplo_net$edges)))
message(sprintf("  country network   : %d nodes, %d links",
                nrow(country_net$nodes), nrow(country_net$links)))

saveRDS(list(haplo = haplo_net, country = country_net, built = Sys.time()),
        file.path(DATA_DIR, "networks.rds"), compress = "xz")

# ================================================================================
# 4. PUBMED SNAPSHOT
# ================================================================================
message("\n[4/4] PubMed snapshot")

# Seeds the Literature tab so the first visitor sees results instantly instead of
# waiting 20-40 s for the live fetch. The app still refreshes in the background
# and the tab's Refresh button still forces a live fetch.
old_wd <- getwd()
setwd("Dashboard")
source("pubmed_module.R")
setwd(old_wd)

papers <- tryCatch(.get_pubmed_papers(force = TRUE), error = function(e) {
  message("  WARNING: PubMed fetch failed (", conditionMessage(e), ")")
  message("  The app will fall back to fetching live on first visit.")
  NULL
})

if (!is.null(papers) && nrow(papers) > 0) {
  message(sprintf("  cached: %d papers", nrow(papers)))
  saveRDS(list(data = papers, fetched = Sys.time()),
          file.path(DATA_DIR, "pubmed_snapshot.rds"), compress = "xz")
} else {
  message("  no snapshot written")
}

# ================================================================================
message("\nDone. Cache files in ", DATA_DIR, ":")
for (f in c("africa_railways.rds", "shipping_routes.rds",
            "networks.rds", "pubmed_snapshot.rds")) {
  p <- file.path(DATA_DIR, f)
  if (file.exists(p)) message(sprintf("  %-24s %s", f, fmt(file.info(p)$size)))
}
