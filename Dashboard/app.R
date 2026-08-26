# ---- Packages ----
library(shiny)
library(bslib)
library(shinycssloaders)
library(leaflet)
library(leaflet.minicharts)
library(tidyverse)
library(tidygeocoder)
library(htmlwidgets)
library(htmltools)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(readxl)
library(visNetwork)
library(ape)
library(stringr)
library(httr)
library(xml2)

local({
  for (pkg in c("igraph",      # visIgraphLayout(); centrality in make_transnet()
                "castor",      # asr_max_parsimony() in make_transnet()
                "hash",        # hash() in make_transnet()
                "plyr",        # count() in make_transnet()
                "tibble",      # rowid_to_column() in make_transnet()
                "magrittr",    # %>% in strainhub_transnet.R
                "tidytree")) { # tidytree::drop.tip()
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Required package '", pkg, "' is not installed.")
    }
  }
})

source("pubmed_module.R")

# make_transnet(), vendored from strainhub (Apache 2.0).
# Strainhub_transnet.R documents the full reasoning.
source("strainhub_transnet.R")

# ---- Configuration ----
MAPBOX_TOKEN      <- "pk.eyJ1Ijoic2F5YWxnIiwiYSI6ImNtaWdha3Y4eDA1YmczZXEybjZvZjE0YTQifQ.T_HqxBxbmcdViI2L0LIzgQ"
MAPBOX_STYLE      <- "https://api.mapbox.com/styles/v1/sayalg/cmkr2q0c3000q01s87umrehau/tiles/256/{z}/{x}/{y}@2x?access_token={accessToken}"
INVALID_LOCATIONS <- c("Laboratory")

# ---- Basemap ----------------------------------------------------------------
# CARTO now requires an API key for the raster basemap tiles at
# basemaps.cartocdn.com, which is what providers$CartoDB.Positron points at.
# Without a key, CARTO serves tiles stamped "requires an API key" — visible as
# soon as you zoom past the levels their CDN still has cached.
#
# A key is free (fair use: 5M tile requests/month) and needs no CARTO account:
#     https://carto.com/basemaps/apikey
#
# Set it here, or leave it empty and set the CARTO_API_KEY environment variable
# (on shinyapps.io: rsconnect::deployApp(..., envVars = "CARTO_API_KEY")).
#
# With a key   -> CartoDB Positron, exactly as before, usable to zoom 20.
# Without one  -> Esri Light Gray Canvas, no key required, visually similar.
#                 Esri only renders that layer to zoom 16, so maxNativeZoom is
#                 pinned there and Leaflet upscales beyond it. Without that pin
#                 Esri returns a "no data here" placeholder tile at z17+, which
#                 would just be a different broken-looking map.
#
# Note CARTO have said the raster basemaps are on a retirement path, so the
# keyless fallback is the more durable option of the two.
CARTO_API_KEY <- Sys.getenv("CARTO_API_KEY", "")

# Esri's World Imagery (the Satellite layer) has its own ceiling: real tiles to
# zoom 19, the same placeholder beyond.
ESRI_IMAGERY_MAX_NATIVE_ZOOM <- 19
ESRI_GRAY_MAX_NATIVE_ZOOM    <- 16
MAP_MAX_ZOOM                 <- 19

# ---- Data paths -------------------------------------------------------------
DATA_DIR <- "data"

SAMPLES_PATH  <- file.path(DATA_DIR, "FinalCOI_metadata.csv")
METADATA_PATH <- file.path(DATA_DIR, "FinalCOI_metadata_w_coor.csv")
INVASIVE_PATH <- file.path(DATA_DIR, "MTM_INVASIVE_VECTOR_SPECIES_20251205.xlsx")
TREE_PATH     <- file.path(DATA_DIR, "COI_subsampled20_per_country_year.aligned.raxml.bestTree")
RAILWAYS_PATH <- file.path(DATA_DIR, "africa_railways.rds")
SHIPPING_PATH <- file.path(DATA_DIR, "shipping_routes.rds")

# Fail fast and loudly at startup rather than with an opaque error mid-render.
local({
  required <- c(SAMPLES_PATH, METADATA_PATH, INVASIVE_PATH, TREE_PATH)
  absent   <- required[!file.exists(required)]
  if (length(absent) > 0) {
    stop("Missing required data file(s): ", paste(absent, collapse = ", "),
         "\nThe app must be run with Dashboard/ as the working directory.")
  }
})

# ================================================================================
# DATA LOADING  (runs once per R process, shared by all sessions)
# ================================================================================

samples_df <- read.csv(SAMPLES_PATH, sep = ",", header = TRUE)
samples_df[samples_df == ""] <- NA

invasive_status <- read_excel(INVASIVE_PATH, sheet = "Data")
invasive_status[invasive_status == ""] <- NA

samples_df <- samples_df %>%
  filter(
    !str_detect(City.Town,    paste(INVALID_LOCATIONS, collapse = "|")) | is.na(City.Town),
    !str_detect(County.State, paste(INVALID_LOCATIONS, collapse = "|")) | is.na(County.State)
  )

# ---- Geocode Missing Coordinates ----
# NOTE: This function makes one API call per missing row and will be extremely
# slow on large datasets. It is kept here for reference only. Run it once
# offline, save the result, then load the pre-geocoded CSV below.
fill_missing_coords <- function(df,
                                country_col = "Country",
                                state_col   = "County.State",
                                city_col    = "City.Town",
                                lat_col     = "Latitude",
                                lon_col     = "Longitude",
                                method      = "osm") {
  df %>%
    mutate(row_id = row_number()) %>%
    group_split(row_id) %>%
    map_dfr(function(row) {
      if (!is.na(row[[lat_col]]) && !is.na(row[[lon_col]])) return(row)
      address <- if (!is.na(row[[city_col]]) && row[[city_col]] != "") {
        str_c(row[[city_col]], row[[state_col]], row[[country_col]], sep = ", ")
      } else {
        row[[country_col]]
      }
      geo <- tryCatch(
        tidygeocoder::geocode(
          tibble(address = address),
          address = address,
          method  = method,
          lat     = "lat",
          long    = "lon",
          quiet   = TRUE
        ),
        error = function(e) NULL
      )
      if (is.null(geo) || is.na(geo$lat) || is.na(geo$lon)) return(row)
      row[[lat_col]] <- geo$lat
      row[[lon_col]] <- geo$lon
      row
    }) %>%
    select(-row_id)
}

# Map data: loaded with base read.csv to match v2 semantics exactly (empty
# strings are kept as "", not coerced to NA).
df_geo   <- read.csv(METADATA_PATH, sep = ",", header = TRUE)
df_haplo <- df_geo %>% filter(!is.na(Haplotype))

# Network data: loaded ONCE with readr::read_csv (empty strings -> NA), matching
# v2's per-render reads. This single object is shared by BOTH network builders
# below, replacing the two separate disk reads v2 performed on every render.
metadata_full <- readr::read_csv(METADATA_PATH, col_names = TRUE, show_col_types = FALSE)

# ---- Cached world polygons (loaded once) ----
world_polygons <- ne_countries(returnclass = "sf") %>%
  mutate(name_clean = tolower(admin))

# ---- Africa railways (cache built by Tools/build_app_cache.R) ----
# One MULTILINESTRING per railway type (rail / light_rail / narrow_gauge /
# preserved), simplified to 500 m and rounded to 4 dp. The raw OSM extract is
# 32,564 separate features and ~16 MB of GeoJSON per browser; this is 3 features
# and ~2.6 MB, which is the difference between a map that pans and one that
# doesn't. Rebuild with: Rscript Tools/build_app_cache.R
africa_railways <- if (file.exists(RAILWAYS_PATH)) {
  readRDS(RAILWAYS_PATH)
} else {
  warning("africa_railways.rds not found — railway layer will be hidden. ",
          "Run: Rscript Tools/build_app_cache.R")
  NULL
}

# ---- Shipping routes (cache built by Tools/build_app_cache.R) ----
# list(routes = sf(route_id, from, to, geometry),
#      route_years = data.frame(route_id, year),
#      years = integer vector)
#
# The raw data repeated each from-to pair once per year — 34,735 features where
# only 1,774 are distinct, all drawn on top of each other. Geometry is now stored
# once and the year index says which routes were active when.
shipping_cache <- if (file.exists(SHIPPING_PATH)) {
  readRDS(SHIPPING_PATH)
} else {
  warning("shipping_routes.rds not found — shipping layer disabled. ",
          "Run: Rscript Tools/build_app_cache.R")
  NULL
}

shipping_routes <- if (!is.null(shipping_cache)) shipping_cache$routes      else NULL
shipping_index  <- if (!is.null(shipping_cache)) shipping_cache$route_years else NULL
shipping_years  <- if (!is.null(shipping_cache)) shipping_cache$years       else integer(0)

# ---- Pre-built networks (cache built by Tools/build_app_cache.R) ----
# Both node/edge tables are computed offline, so no session pays for the
# parsimony reconstruction or the centrality maths.
NETWORKS_PATH <- file.path(DATA_DIR, "networks.rds")
networks_cache <- if (file.exists(NETWORKS_PATH)) {
  readRDS(NETWORKS_PATH)
} else {
  warning("networks.rds not found — networks will be computed at runtime. ",
          "Run: Rscript Tools/build_app_cache.R")
  NULL
}

# ---- PubMed snapshot (cache built by Tools/build_app_cache.R) ----
# Seeds the Literature tab so the first visitor is not left waiting 20-40 s for
# the live NCBI fetch. Optional — the module falls back to fetching live.
PUBMED_SNAPSHOT_PATH <- file.path(DATA_DIR, "pubmed_snapshot.rds")
if (!seed_pubmed_cache(PUBMED_SNAPSHOT_PATH)) {
  message("No PubMed snapshot found — the Literature tab will fetch live on ",
          "first visit. Run: Rscript Tools/build_app_cache.R")
}

# ================================================================================
# COLOR UTILITIES
# ================================================================================

create_high_contrast_palette <- function(n) {
  base_colors <- c(
    "#E41A1C", "#4DAF4A", "#377EB8", "#FF7F00", "#984EA3",
    "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
    "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
    "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666"
  )
  if (n <= length(base_colors)) return(base_colors[seq_len(n)])
  colorRampPalette(base_colors)(n)
}

create_haplotype_palette <- function(data, category_col = "Haplotype") {
  all_haplotypes <- data %>%
    filter(!is.na(.data[[category_col]])) %>%
    pull(.data[[category_col]]) %>%
    unique()
  numeric_parts <- as.numeric(str_extract(all_haplotypes, "\\d+"))
  if (!all(is.na(numeric_parts))) {
    all_haplotypes <- all_haplotypes[order(numeric_parts, na.last = TRUE)]
  }
  colors <- create_high_contrast_palette(length(all_haplotypes))
  setNames(as.character(colors), all_haplotypes)
}

# Palette depends only on static data, so compute it once at startup rather than
# once per session.
haplotype_palette <- create_haplotype_palette(df_haplo, category_col = "Haplotype")

# ================================================================================
# MAP FUNCTIONS
# ================================================================================

# The light "Base" layer, keyed or keyless depending on CARTO_API_KEY above.
add_base_tiles <- function(map) {
  if (nzchar(CARTO_API_KEY)) {
    addTiles(
      map,
      urlTemplate = paste0(
        "https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}.png?key=",
        CARTO_API_KEY
      ),
      attribution = paste(
        '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a>',
        'contributors &copy; <a href="https://carto.com/attributions">CARTO</a>'
      ),
      options = tileOptions(subdomains = "abcd", maxZoom = 20),
      group   = "Base"
    )
  } else {
    addProviderTiles(
      map,
      providers$Esri.WorldGrayCanvas,
      group   = "Base",
      options = providerTileOptions(maxNativeZoom = ESRI_GRAY_MAX_NATIVE_ZOOM,
                                    maxZoom       = MAP_MAX_ZOOM)
    )
  }
}

create_base_map <- function(tile_size = 256) {
  leaflet(options = leafletOptions(maxZoom = MAP_MAX_ZOOM)) %>%
    addMapPane("countryStatus", zIndex = 200) %>%
    addMapPane("charts",        zIndex = 400) %>%
    add_base_tiles() %>%
    addProviderTiles(
      providers$Esri.WorldImagery,
      group   = "Satellite",
      options = providerTileOptions(maxNativeZoom = ESRI_IMAGERY_MAX_NATIVE_ZOOM,
                                    maxZoom       = MAP_MAX_ZOOM)
    ) %>%
    addTiles(
      urlTemplate = MAPBOX_STYLE,
      options     = tileOptions(accessToken = MAPBOX_TOKEN, tileSize = tile_size),
      group       = "Land Use"
    ) %>%
    addLayersControl(
      baseGroups    = c("Base", "Satellite"),
      overlayGroups = c("Land Use", "Country Status", "Railways", "Shipping Routes"),
      options       = layersControlOptions(collapsed = FALSE)
    ) %>%
    hideGroup("Railways") %>%
    hideGroup("Shipping Routes") %>%
    onRender("
      function(el, x) {
        var map = this;
        var attached = {};
        // Poll for freshly-drawn minicharts and wire a click handler to each.
        // Charts are re-drawn whenever the temporal filter changes, so the
        // interval keeps running, but it only touches charts it hasn't seen.
        map._minichartClickPoll = setInterval(function() {
          if (!map._minicharts) return;
          map._minicharts.forEach(function(chart) {
            var layerId = chart.options.layerId;
            if (!layerId || attached[layerId]) return;
            var chartEl = chart._container;
            if (!chartEl) return;
            attached[layerId] = true;
            chartEl.addEventListener('click', function(e) {
              Shiny.setInputValue('map_marker_click', {id: layerId}, {priority: 'event'});
            }, false);
          });
        }, 500);
      }
    ")
}

add_pie_markers <- function(map, data,
                            lat_col       = "Latitude",
                            lon_col       = "Longitude",
                            category_col  = "Haplotype",
                            city_col      = "City.Town",
                            state_col     = "County.State",
                            country_col   = "Country",
                            radius        = 20,
                            add_legend    = TRUE,
                            color_palette = NULL) {

  city_lookup <- data %>%
    filter(!is.na(.data[[lat_col]]), !is.na(.data[[lon_col]])) %>%
    group_by(.data[[lat_col]], .data[[lon_col]]) %>%
    summarise(
      city = {
        cities <- .data[[city_col]][!is.na(.data[[city_col]]) & .data[[city_col]] != ""]
        if (length(cities) > 0) {
          names(which.max(table(cities)))
        } else {
          states <- .data[[state_col]][!is.na(.data[[state_col]]) & .data[[state_col]] != ""]
          if (length(states) > 0) {
            names(which.max(table(states)))
          } else {
            countries <- .data[[country_col]][!is.na(.data[[country_col]]) & .data[[country_col]] != ""]
            if (length(countries) > 0) names(which.max(table(countries))) else "Unknown"
          }
        }
      },
      .groups = "drop"
    )

  df_pie <- data %>%
    filter(!is.na(.data[[lat_col]]), !is.na(.data[[lon_col]])) %>%
    count(.data[[lat_col]], .data[[lon_col]], .data[[category_col]]) %>%
    pivot_wider(names_from = all_of(category_col), values_from = n, values_fill = 0) %>%
    left_join(city_lookup, by = c(lat_col, lon_col))

  coords <- df_pie %>% select(all_of(c(lat_col, lon_col)))
  pies   <- df_pie %>% select(-all_of(c(lat_col, lon_col, "city")))

  numeric_parts <- as.numeric(str_extract(colnames(pies), "\\d+"))
  if (!all(is.na(numeric_parts))) pies <- pies[, order(numeric_parts, na.last = TRUE)]

  if (!is.null(color_palette)) {
    colors <- unname(color_palette[colnames(pies)])
    if (any(is.na(colors))) {
      warning("Some haplotypes not found in color_palette; using fallback colors.")
      colors <- create_high_contrast_palette(ncol(pies))
    }
  } else {
    colors <- create_high_contrast_palette(ncol(pies))
  }
  colors <- as.character(colors)

  layer_ids <- make.unique(df_pie$city)
  color_map <- setNames(colors, colnames(pies))
  popup_list <- setNames(vector("list", nrow(df_pie)), layer_ids)

  for (i in seq_len(nrow(df_pie))) {
    city_name <- htmltools::htmlEscape(df_pie$city[i])
    hap_vals  <- as.numeric(pies[i, ])
    names(hap_vals) <- colnames(pies)
    present <- sort(hap_vals[hap_vals > 0], decreasing = TRUE)
    total   <- sum(present)

    rows_html <- paste0(
      "<tr>",
      "<td style='padding:2px 8px 2px 0'>",
      "<span style='display:inline-block;width:10px;height:10px;background-color:",
      color_map[names(present)], ";border-radius:50%;margin-right:5px'></span>",
      names(present), "</td>",
      "<td style='padding:2px 0'>", present,
      " <span style='color:#888'>(", round(100 * present / total, 1), "%)</span></td>",
      "</tr>",
      collapse = ""
    )

    popup_list[[i]] <- paste0(
      "<div style='font-family:sans-serif;min-width:160px'>",
      "<b style='font-size:13px'>", city_name, "</b>",
      "<div style='color:#666;font-size:11px;margin:3px 0 6px'>",
      total, " sample", if (total != 1) "s" else "", "</div>",
      "<table style='font-size:12px;border-collapse:collapse'>",
      rows_html, "</table></div>"
    )
  }

  popup_html <- unlist(popup_list[layer_ids])

  map_with_charts <- addMinicharts(
    map          = map,
    lng          = coords[[lon_col]],
    lat          = coords[[lat_col]],
    type         = "pie",
    chartdata    = as.matrix(pies),
    width        = radius,
    height       = radius,
    colorPalette = colors,
    layerId      = layer_ids,
    popup        = popupArgs(html = popup_html),
    showLabels   = FALSE,
    legend       = FALSE
  )

  if (add_legend) {
    map_with_charts <- map_with_charts %>%
      addLegend(
        position = "bottomleft",
        colors   = colors,
        labels   = colnames(pies),
        title    = "Haplotypes",
        opacity  = 1,
        layerId  = "haplotype_legend"
      )
  }

  # Add invisible markers on top of pie charts to capture clicks
  map_with_charts <- map_with_charts %>%
    addMarkers(
      lng      = coords[[lon_col]],
      lat      = coords[[lat_col]],
      layerId  = layer_ids,
      popup    = unname(popup_html),
      icon     = makeIcon(
        iconUrl = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNk+M9QDwADhgGAWjR9awAAAABJRU5ErkJggg==",
        iconWidth = 60,
        iconHeight = 60
      ),
      options = markerOptions(opacity = 0, zIndexOffset = 1000)
    )

  map_with_charts
}

add_country_status_layer <- function(map, data,
                                     country_col = "Country",
                                     status_col  = "Native.Invasive",
                                     world       = world_polygons) {
  country_status <- data %>%
    mutate(
      country_clean = tolower(.data[[country_col]]),
      status_clean  = tolower(.data[[status_col]])
    ) %>%
    filter(!is.na(country_clean), !is.na(status_clean)) %>%
    count(country_clean, status_clean) %>%
    slice_max(n, n = 1, by = country_clean) %>%
    select(country_clean, status_clean)

  world_status <- world %>%
    left_join(country_status, by = c("name_clean" = "country_clean"))

  pal <- colorFactor(
    palette  = c("native" = "#4CAF50", "invasive" = "#FF9800"),
    domain   = c("native", "invasive"),
    na.color = "transparent"
  )

  map %>%
    addPolygons(
      data        = world_status,
      fillColor   = ~pal(status_clean),
      fillOpacity = 0.15,
      color       = "#444444",
      weight      = 0.5,
      opacity     = 0.5,
      group       = "Country Status",
      options     = pathOptions(pane = "countryStatus"),
      label       = ~paste0(admin, ": ", if_else(is.na(status_clean), "No Data", status_clean))
    ) %>%
    addLegend(
      position = "bottomright",
      pal      = pal,
      values   = c("native", "invasive"),
      title    = "Native vs. Invasive",
      group    = "Country Status"
    )
}

add_chart_interactions <- function(map) {
  map %>%
    onRender("
      function(el, x) {
        var map = this;
        if (!map._minicharts) return;
        map._minicharts.forEach(function(chart) {
          var chartEl = chart._container;
          if (!chartEl) return;
          chartEl.style.position = 'absolute';
          chartEl.style.zIndex   = 400;
          map.getPane('charts').appendChild(chartEl);
        });
      }
    ")
}

add_railway_layer <- function(map, railways = africa_railways) {
  if (is.null(railways) || nrow(railways) == 0) return(map)

  # Color by railway type
  pal <- colorFactor(
    palette = c(
      "rail"          = "#4a4a4a",
      "light_rail"    = "#8a6db5",
      "narrow_gauge"  = "#c0622a",
      "preserved"     = "#2e7d52"
    ),
    domain  = c("rail", "light_rail", "narrow_gauge", "preserved"),
    na.color = "#888888"
  )

  map %>%
    addMapPane("railways", zIndex = 300) %>%
    addPolylines(
      data        = railways,
      color       = ~pal(fclass),
      weight      = 1.2,
      opacity     = 0.75,
      group       = "Railways",
      options     = pathOptions(pane = "railways"),
      # The cache merges every railway of a type into one MULTILINESTRING, so
      # there is no per-railway name to show any more. That trade buys a ~6x
      # smaller payload and takes Leaflet from 32,564 SVG paths down to 3, which
      # is what was making the map crawl. To get individual names back, drop the
      # group_by(fclass)/summarise() step in Tools/build_app_cache.R and keep the
      # `name` column — expect the slowdown to come with it.
      label        = ~paste0("Railway type: ", fclass),
      labelOptions = labelOptions(style = list("font-size" = "12px"))
    ) %>%
    addLegend(
      position = "bottomright",
      pal      = pal,
      values   = railways$fclass,
      title    = "Railway type",
      group    = "Railways",
      layerId  = "railway_legend"
    )
}

add_shipping_layer <- function(map, routes_sf, highlight_location = NULL) {
  map <- clearGroup(map, "Shipping Routes")
  if (is.null(routes_sf) || nrow(routes_sf) == 0) return(map)

  # Add the map pane
  map <- map %>% addMapPane("shipping", zIndex = 250)

  if (!is.null(highlight_location) && highlight_location != "") {
    # Split routes into highlighted and non-highlighted
    routes_sf$is_highlighted <- (routes_sf$from == highlight_location |
                                   routes_sf$to == highlight_location)

    # Add non-highlighted routes first (low opacity)
    non_highlighted <- routes_sf[!routes_sf$is_highlighted, ]
    if (nrow(non_highlighted) > 0) {
      map <- map %>%
        addPolylines(
          data    = non_highlighted,
          color   = "#0077b6",
          weight  = 1,
          opacity = 0.1,
          group   = "Shipping Routes",
          options = pathOptions(pane = "shipping")
        )
    }

    # Add highlighted routes on top (high opacity)
    highlighted <- routes_sf[routes_sf$is_highlighted, ]
    if (nrow(highlighted) > 0) {
      map <- map %>%
        addPolylines(
          data    = highlighted,
          color   = "#d62728",
          weight  = 2.5,
          opacity = 0.85,
          group   = "Shipping Routes",
          options = pathOptions(pane = "shipping")
        )
    }
  } else {
    # No selection: show all routes with default styling
    map <- map %>%
      addPolylines(
        data    = routes_sf,
        color   = "#0077b6",
        weight  = 1,
        opacity = 0.4,
        group   = "Shipping Routes",
        options = pathOptions(pane = "shipping")
      )
  }

  map
}

# ================================================================================
# NETWORK BUILDERS  (served from the offline cache; computed only as a fallback)
# ================================================================================
# Neither network depends on user input, so both are built offline by
# Tools/build_app_cache.R and read straight out of networks.rds. That keeps the
# parsimony reconstruction, the centrality maths and the tree parsing off the
# server entirely — and means castor/hash/plyr/igraph are only touched if the
# cache is missing.
#
# The runtime builders below are kept as a fallback so the app still works on a
# checkout where the cache has not been generated yet.

.app_cache <- new.env(parent = emptyenv())

# ---- Panel 1: Strainhub country transmission network ----
build_haplo_network_data <- function() {
  if (!is.null(.app_cache$haplo_net)) return(.app_cache$haplo_net)

  if (!is.null(networks_cache$haplo)) {
    .app_cache$haplo_net <- networks_cache$haplo
    return(.app_cache$haplo_net)
  }

  result <- tryCatch({
    treedata <- ape::read.tree(TREE_PATH)
    metadata <- metadata_full   # reuse the already-loaded metadata

    ## Adjust column names
    names(metadata)[names(metadata) == "FASTA.ID"]     <- "Accession"
    names(metadata)[names(metadata) == "Accession.ID"] <- "Unique"

    ## Prune Tree
    taxa_without_data <- setdiff(treedata$tip.label, metadata$Accession)
    pruned_tree <- tidytree::drop.tip(treedata, taxa_without_data)

    ## Create Country with Haplotype column
    metadata$haplo_concat <- paste(metadata$Country, metadata$Haplotype, sep = "_")

    ## Make the Transmission Network
    # metricsOutputFile = "" is required for deployment: the default
    # ("StrainHub_metrics.csv") makes make_transnet() write into the working
    # directory, which is read-only on shinyapps.io / Posit Connect and would
    # abort the render. "" tells strainhub to skip the file entirely.
    graph <- make_transnet(pruned_tree,
                           metadata,
                           columnSelection   = "Country",
                           centralityMetric  = 3,   # Betweenness Centrality
                           treeType          = "parsimonious",
                           metricsOutputFile = "")

    nodes <- graph$x$nodes %>%
      mutate(shape = "dot",
             color = "grey",
             font.size = 35)

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
  }, error = function(e) {
    warning("build_haplo_network_data() failed: ", conditionMessage(e))
    NULL
  })

  .app_cache$haplo_net <- result
  result
}

# ---- Panel 2: Country connectivity (shared haplotypes) network ----
build_country_network_data <- function() {
  if (!is.null(.app_cache$country_net)) return(.app_cache$country_net)

  if (!is.null(networks_cache$country)) {
    .app_cache$country_net <- networks_cache$country
    return(.app_cache$country_net)
  }

  result <- tryCatch({
    # Clean metadata (reuse the already-loaded object)
    metadata_clean_cn <- metadata_full %>%
      filter(!is.na(Country) & Country != "" & Country != "NA") %>%
      filter(!is.na(Haplotype) & Haplotype != "" & Haplotype != "NA")

    # Random subsampling (20 per Country/Year)
    set.seed(123)
    metadata_subsampled_cn <- metadata_clean_cn %>%
      group_by(Country, Year) %>%
      slice_sample(n = 20) %>%
      ungroup()

    # Binary haplotype matrix
    h_matrix_raw    <- table(metadata_subsampled_cn$Haplotype, metadata_subsampled_cn$Country)
    h_matrix_binary <- h_matrix_raw
    h_matrix_binary[h_matrix_binary > 0] <- 1

    # Shared unique haplotype count between countries
    country_matrix <- as.matrix(t(h_matrix_binary) %*% h_matrix_binary)
    diag(country_matrix) <- 0

    # Nodes
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

    # Edges (upper triangle only)
    edges_idx <- which(country_matrix > 0, arr.ind = TRUE)
    edges_idx <- edges_idx[edges_idx[, 1] < edges_idx[, 2], , drop = FALSE]
    links_country <- data.frame(
      from  = edges_idx[, 1],
      to    = edges_idx[, 2],
      value = country_matrix[edges_idx]
    )

    list(nodes = nodes_country, links = links_country)
  }, error = function(e) {
    warning("build_country_network_data() failed: ", conditionMessage(e))
    NULL
  })

  .app_cache$country_net <- result
  result
}

# ================================================================================
# THEME
# ================================================================================

app_theme <- bs_theme(
  version       = 5,
  bg            = "#f4f6f9",
  fg            = "#1c2b3a",
  primary       = "#1a4f72",       # deep ocean blue  — navbar, buttons
  secondary     = "#2e86ab",       # mid blue         — accents
  success       = "#4CAF50",       # native green     — matches map
  warning       = "#FF7F00",       # invasive orange  — matches map
  base_font     = font_google("Inter"),
  heading_font  = font_google("Merriweather"),
  code_font     = font_google("Fira Code"),
  "navbar-bg"              = "#1a4f72",
  "navbar-brand-color"     = "#ffffff",
  "navbar-light-color"     = "#ffffff",
  "navbar-dark-color"      = "#ffffff",
  "nav-link-color"         = "rgba(255,255,255,0.85)",
  "nav-link-hover-color"   = "#ffffff",
  "card-border-radius"     = "0.6rem",
  "card-box-shadow"        = "0 2px 8px rgba(0,0,0,0.08)"
)

# ================================================================================
# UI
# ================================================================================

ui <- page_navbar(
  title = tags$span(
    tags$img(src = "logo.png", height = "100px", style = "margin-right:10px;vertical-align:middle"),
    "The Vector Invasion Atlas"
  ),
  theme    = app_theme,
  id       = "main_nav",
  fillable = FALSE,

  # Link the external stylesheet
  header = tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
    tags$link(rel = "icon", type = "image/png", href = "favicon.png")
  ),

  # ---- Gene selector (top-right of navbar) -----------------------------------
  nav_item(
    tags$div(
      style = "padding: 6px 10px;",
      selectInput(
        "selected_gene",
        label    = NULL,
        choices  = c("COI" = "COI"),   # add more genes here later
        selected = "COI",
        width    = "120px"
      )
    )
  ),

  # ---- Tab 1: Interactive Map ------------------------------------------------
  nav_panel(
    title = tags$span(icon("globe"), "Distribution Map"),
    value = "map_tab",

    layout_sidebar(
      fillable = FALSE,

      sidebar = sidebar(
        width = 280,
        open  = TRUE,

        tags$div(class = "sidebar-section-header", "Temporal Filter"),

        radioButtons(
          "time_mode", label = NULL,
          choices  = c("Exact Year"   = "exact",
                       "Up to Year"   = "cumulative",
                       "Year Range"   = "range",
                       "All Years"    = "all"),
          selected = "all"
        ),

        conditionalPanel(
          condition = "input.time_mode === 'exact' || input.time_mode === 'cumulative'",
          sliderInput(
            "year_filter", "Year:",
            min     = 1990,
            max     = 2024,
            value   = 2024,
            step    = 1,
            sep     = "",
            animate = animationOptions(interval = 1000, loop = FALSE)
          )
        ),

        conditionalPanel(
          condition = "input.time_mode === 'range'",
          sliderInput(
            "year_range", "Year range:",
            min   = 1990,
            max   = 2024,
            value = c(2016, 2024),
            step  = 1,
            sep   = ""
          )
        ),

        hr(class = "sidebar-divider"),
        tags$div(class = "sidebar-section-header", "Summary"),
        uiOutput("sample_count_ui"),
        uiOutput("no_data_warning")
      ),

      # Map fills available height; spinner while loading
      card(
        full_screen = TRUE,
        class       = "map-card",
        withSpinner(
          leafletOutput("map", height = "780px"),
          color = "#1a4f72"
        )
      )
    )
  ),

  # ---- Tab 2: Phylogenetic Tree -----------------------------------------------
  nav_panel(
    title = tags$span(icon("tree"), "Phylogenetic Tree"),
    value = "tree_tab",
    layout_columns(
      col_widths = c(3, 9),

      # Left: info card (update description as needed)
      card(
        card_header(tags$span(icon("circle-info"), " About this Tree")),
        class = "controls-card",
        tags$p(class = "help-text",
               "Maximum likelihood phylogenetic tree of", tags$em("Anopheles stephensi"),
               "COI haplotypes inferred with RAxML. Tree is hosted and visualised
              interactively via iTOL. Use the iTOL toolbar to zoom, collapse clades,
              and export the figure.")
      ),

      # Right: iTOL interactive viewer
      card(
        full_screen = TRUE,
        card_header(
          "iTOL — COI Maximum Likelihood Tree",
          tags$a(
            href     = "tree.svg",
            download = "anopheles_stephensi_tree.svg",
            class    = "btn btn-sm btn-outline-primary float-end",
            icon("download"), " Download SVG"
          )
        ),
        tags$img(src = "tree.svg", width = "100%", style = "height:700px; object-fit:contain;")
      )
    )
  ),

  # ---- Tab 3: Gene Flow Network -----------------------------------------------
  nav_panel(
    title = tags$span(icon("circle-nodes"), "Gene Network"),
    value = "network_tab",

    layout_columns(
      col_widths = c(3, 9),

      # Left: controls
      card(
        hr(class = "sidebar-divider"),
        tags$p(class = "help-text",
               style = "padding: 72px; font-size: 16px; margin: 8px 0;",
               icon("circle-info"),
               " Node size reflects haplotype frequency. Edge direction reflects inferred character state changes.")
      ),

      # Right: network
      card(
        full_screen = TRUE,
        card_header("Strainhub: Haplotype Gene Flow Network"),
        withSpinner(
          visNetworkOutput("haplo_network", height = "700px"),
          color = "#1a4f72"
        )
      )
    ),

    layout_columns(
      col_widths = c(3, 9),

      # Left: controls for country connectivity
      card(
        hr(class = "sidebar-divider"),
        tags$p(class = "help-text",
               style = "padding: 72px; font-size: 16px; margin: 8px 0;",
               icon("circle-info"),
               " Node size reflects number of included samples (capped at 20 per country/year). Edge thickness reflects the number of shared unique haplotypes between countries.")
      ),

      # Right: country connectivity network
      card(
        full_screen = TRUE,
        card_header("An. stephensi: Country Connectivity"),
        withSpinner(
          visNetworkOutput("country_network", height = "700px"),
          color = "#1a4f72"
        )
      )
    )
  ),

  # ---- Tab 4: Literature (PubMed module) --------------------------------------
  nav_panel(
    title = tags$span(icon("book-open"), "Literature"),
    value = "lit_tab",
    pubBrowserUI("pubs")
  ),

  # ---- Footer / About (nav_menu) ----------------------------------------------
  nav_menu(
    title = tags$span(icon("circle-info"), "About"),
    align = "right",

    nav_panel(
      "About this Tool",
      card(
        class = "about-card",
        card_body(
          tags$h4("About VIA"),
          tags$p(
            "This dashboard visualizes COI haplotype diversity and geographic distribution of ",
            tags$em("Anopheles stephensi"), " and related vectors across native and invasive ranges."
          ),
          tags$h5("Data"),
          tags$p(
            "Scripts available here: ", tags$a(href= "https://github.com/sayalg/VIA", "https://github.com/sayalg/VIA")
          ),
          tags$h5("Contact"),
          # TODO: update with your actual contact info
          tags$p("For questions or data requests, contact: ",
                 tags$a(href = "mailto:Tamar_Carter@baylor.edu", "Tamar_Carter@baylor.edu")),
          tags$h5("Citation"),
          tags$p(tags$em("Coming Soon")),
          tags$hr(),
          tags$p(class = "text-muted small",
                 paste0("Built with R Shiny & bslib | Last updated: ", format(Sys.Date(), "%B %Y")))
        )
      )
    )
  )
)

# ================================================================================
# SERVER
# ================================================================================

server <- function(input, output, session) {

  # ---- Update slider range from actual data ----
  available_years <- df_haplo %>%
    filter(!is.na(Year)) %>%
    pull(Year) %>%
    unique() %>%
    sort()

  updateSliderInput(session, "year_filter",
                    min   = min(available_years),
                    max   = max(available_years),
                    value = max(available_years)
  )
  updateSliderInput(session, "year_range",
                    min   = min(available_years),
                    max   = max(available_years),
                    value = c(min(available_years), max(available_years))
  )

  # Track selected location for shipping route highlighting
  selected_location <- reactiveVal(NULL)

  # Create city-to-country mapping
  city_to_country <- df_haplo %>%
    filter(!is.na(City.Town), !is.na(Country)) %>%
    select(City.Town, Country) %>%
    distinct() %>%
    deframe()

  # ==============================================================================
  # TAB 1 — MAP
  # ==============================================================================

  filtered_data <- reactive({
    switch(input$time_mode,
           "exact"      = df_haplo %>% filter(Year == input$year_filter),
           "cumulative" = df_haplo %>% filter(Year <= input$year_filter),
           "range"      = df_haplo %>% filter(Year >= input$year_range[1],
                                              Year <= input$year_range[2]),
           "all"        = df_haplo
    )
  })

  output$sample_count_ui <- renderUI({
    n <- nrow(filtered_data())
    s <- if (n != 1) "s" else ""
    label <- switch(input$time_mode,
                    "exact"      = paste0(n, " sample", s, " from ", input$year_filter),
                    "cumulative" = paste0(n, " sample", s, " up to ", input$year_filter),
                    "range"      = paste0(n, " sample", s, ", ", input$year_range[1],
                                          "–", input$year_range[2]),
                    "all"        = paste0(n, " sample", s, " (all years)")
    )
    tags$div(class = "sample-count-badge", label)
  })

  output$no_data_warning <- renderUI({
    if (nrow(filtered_data()) == 0 && input$time_mode != "all") {
      year_label <- if (input$time_mode == "range") {
        paste(input$year_range[1], "–", input$year_range[2])
      } else {
        input$year_filter
      }
      tags$div(
        class = "alert alert-warning mt-2 p-2 small",
        icon("triangle-exclamation"),
        paste(" No samples found for", year_label)
      )
    }
  })

  # Base map with static haplotype legend — renders once.
  #
  # NOTE: the Railways layer is deliberately NOT added here. It starts hidden
  # (hideGroup in create_base_map), and adding it at render time meant every
  # visitor downloaded the full railway geometry before they could see anything,
  # whether or not they ever switched the layer on. It is now drawn on demand by
  # the observer below.
  output$map <- renderLeaflet({
    create_base_map() %>%
      add_country_status_layer(
        data        = invasive_status,
        country_col = "COUNTRY_NAME",
        status_col  = "INVASIVE_STATUS"
      ) %>%
      addLegend(
        position = "bottomleft",
        colors   = haplotype_palette,
        labels   = names(haplotype_palette),
        title    = "Haplotypes",
        opacity  = 1,
        layerId  = "haplotype_legend"
      )
  })

  # ---- Lazy heavy layers ------------------------------------------------------
  # input$map_groups holds the overlay groups the user currently has ticked.
  # Both heavy layers are hidden by default, so gating on this means a first page
  # load ships none of their geometry.
  railways_drawn <- reactiveVal(FALSE)

  layer_on <- function(name) {
    isTRUE(name %in% input$map_groups)
  }

  observeEvent(input$map_groups, {
    if (layer_on("Railways") && !railways_drawn()) {
      leafletProxy("map") %>% add_railway_layer()
      railways_drawn(TRUE)   # geometry is static; draw it at most once per session
    }
  }, ignoreNULL = TRUE)

  # Update pie charts on filter change only
  observe({
    req(filtered_data())
    fd    <- filtered_data()
    proxy <- leafletProxy("map") %>% clearMinicharts()
    if (nrow(fd) > 0) {
      proxy <- add_pie_markers(
        map           = proxy,
        data          = fd,
        lat_col       = "Latitude",
        lon_col       = "Longitude",
        category_col  = "Haplotype",
        city_col      = "City.Town",
        state_col     = "County.State",
        country_col   = "Country",
        radius        = 30,
        add_legend    = FALSE,
        color_palette = haplotype_palette
      )
    }
    add_chart_interactions(proxy)
  })

  # Handle pie chart clicks via invisible markers
  observeEvent(input$map_marker_click, {
    click <- input$map_marker_click
    if (is.null(click) || is.null(click$id)) {
      selected_location(NULL)
    } else {
      clicked_city <- click$id
      # Convert city name to country
      country <- city_to_country[clicked_city]
      if (!is.na(country)) {
        selected_location(country)
      }
    }
  })

  # Reset on empty map click
  observeEvent(input$map_click, {
    if (is.null(input$map_click$id)) {
      selected_location(NULL)
    }
  }, ignoreInit = TRUE)

  # Routes matching the current temporal filter. Geometry is stored once per
  # from-to pair and selected through the year index, rather than rbind-ing a
  # separate full copy of every year's geometry together.
  filtered_shipping_sf <- reactive({
    if (is.null(shipping_routes) || length(shipping_years) == 0) return(NULL)

    target_years <- switch(input$time_mode,
                           "exact"      = shipping_years[shipping_years == input$year_filter],
                           "cumulative" = shipping_years[shipping_years <= input$year_filter],
                           "range"      = shipping_years[shipping_years >= input$year_range[1] &
                                                           shipping_years <= input$year_range[2]],
                           "all"        = shipping_years
    )
    if (length(target_years) == 0) return(NULL)

    active <- unique(shipping_index$route_id[shipping_index$year %in% target_years])
    if (length(active) == 0) return(NULL)

    shipping_routes[shipping_routes$route_id %in% active, ]
  })

  # Redraw base routes when the filter changes — but only while the user
  # actually has the layer switched on. This is the single biggest saving in the
  # app: the layer is hidden by default, and previously every session pushed the
  # full route geometry to the browser regardless.
  observe({
    map <- leafletProxy("map") %>%
      clearGroup("Shipping Routes") %>%
      clearGroup("Shipping Highlight")

    if (!layer_on("Shipping Routes")) return()

    routes <- filtered_shipping_sf()
    if (is.null(routes) || nrow(routes) == 0) return()

    map %>%
      addMapPane("shipping", zIndex = 250) %>%
      addPolylines(
        data    = routes,
        color   = "#0077b6",
        weight  = 1,
        opacity = 0.4,
        group   = "Shipping Routes",
        options = pathOptions(pane = "shipping")
      )
  })

  # Update highlight layer only when selection changes
  observe({
    leafletProxy("map") %>% clearGroup("Shipping Highlight")

    if (!layer_on("Shipping Routes")) return()

    routes   <- filtered_shipping_sf()
    location <- selected_location()
    if (is.null(routes) || is.null(location) || location == "") return()
    highlighted <- routes[routes$from == location | routes$to == location, ]
    if (nrow(highlighted) == 0) return()
    leafletProxy("map") %>%
      addPolylines(
        data    = highlighted,
        color   = "#d62728",
        weight  = 2.5,
        opacity = 0.85,
        group   = "Shipping Highlight",
        options = pathOptions(pane = "shipping")
      )
  })

  # ==============================================================================
  # TAB 2 — PHYLOGENETIC TREE
  # ==============================================================================

    # Tree rendered from iTOL in utility section

  # ==============================================================================
  # TAB 3, PANEL 1 — STRAINHUB NETWORK
  # ==============================================================================

  output$haplo_network <- renderVisNetwork({
    net <- build_haplo_network_data()
    validate(need(!is.null(net), "Transmission network could not be built."))

    visNetwork(net$nodes, net$edges,
               main = "Strainhub: Country Transmission") %>%
      visNodes(shape = "ellipse") %>%
      visLayout(randomSeed = 12) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = list(enabled = TRUE)) %>%
      visPhysics(solver = "repulsion", repulsion = list(nodeDistance = 200))
  })

  # ==============================================================================
  # TAB 3, PANEL 2 — HAPLOTYPE NETWORK
  # ==============================================================================

  output$country_network <- renderVisNetwork({
    net <- build_country_network_data()
    validate(need(!is.null(net), "Country connectivity network could not be built."))

    visNetwork(net$nodes, net$links,
               main = "An. stephensi: Country Connectivity") %>%
      visLayout(randomSeed = 1) %>%
      visIgraphLayout(layout = "layout_nicely") %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visEdges(color = "rgba(100, 100, 100, 0.4)", smooth = TRUE) %>%
      visGroups(groupname = "Afghanistan",          color = "#CD9B1D") %>%
      visGroups(groupname = "Djibouti",             color = "#00008B") %>%
      visGroups(groupname = "Ethiopia",             color = "#009ACD") %>%
      visGroups(groupname = "India",                color = "#EE9A00") %>%
      visGroups(groupname = "Iran",                 color = "#CD853F") %>%
      visGroups(groupname = "Kenya",                color = "#6CA6CD") %>%
      visGroups(groupname = "Niger",                color = "#7CCD7C") %>%
      visGroups(groupname = "Pakistan",             color = "#EEDC82") %>%
      visGroups(groupname = "SaudiArabia",          color = "#556B2F") %>%
      visGroups(groupname = "Somalia",              color = "#BCD2EE") %>%
      visGroups(groupname = "Sri Lanka",            color = "#BFEFFF") %>%
      visGroups(groupname = "Sudan",                color = "#1874CD") %>%
      visGroups(groupname = "United Arab Emirates", color = "#FF7F00") %>%
      visGroups(groupname = "Yemen",                color = "#4169E1") %>%
      visPhysics(solver = "repulsion", repulsion = list(nodeDistance = 300))
  })

  # ==============================================================================
  # TAB 4 — LITERATURE (PubMed module)
  # ==============================================================================

  pubBrowserServer("pubs")

}

# ================================================================================
# RUN
# ================================================================================

shinyApp(ui = ui, server = server)
