# ---- Packages ----
# Only what the running app needs. Everything used to *build* the data — geometry
# simplification, the invasive-status spreadsheet, country polygons, StrainHub,
# igraph layouts, geocoding — lives in Tools/build_app_cache.R, so the server
# installs and loads far less and every R process starts faster.
library(shiny)
library(bslib)
library(shinycssloaders)
library(leaflet)
library(leaflet.minicharts)
library(dplyr)
library(tidyr)
library(sf)           # country-status polygons
library(visNetwork)

source("pubmed_module.R")   # also loads stringr, httr, xml2

# ---- Configuration ----
MAPBOX_TOKEN <- "pk.eyJ1Ijoic2F5YWxnIiwiYSI6ImNtaWdha3Y4eDA1YmczZXEybjZvZjE0YTQifQ.T_HqxBxbmcdViI2L0LIzgQ"
MAPBOX_STYLE <- "https://api.mapbox.com/styles/v1/sayalg/cmkr2q0c3000q01s87umrehau/tiles/256/{z}/{x}/{y}@2x?access_token={accessToken}"

# ---- Basemap ----------------------------------------------------------------
# CARTO now requires an API key for basemaps.cartocdn.com raster tiles (free:
# https://carto.com/basemaps/apikey). Without one, CARTO serves tiles stamped
# "requires an API key". Set CARTO_API_KEY (on shinyapps.io:
# rsconnect::deployApp(..., envVars = "CARTO_API_KEY")) for CartoDB Positron;
# otherwise the keyless Esri Light Gray Canvas is used.
#
# Esri only renders Gray Canvas to zoom 16 and World Imagery to 19, returning a
# "no data" placeholder beyond, so maxNativeZoom is pinned and Leaflet upscales.
CARTO_API_KEY <- Sys.getenv("CARTO_API_KEY", "")
ESRI_IMAGERY_MAX_NATIVE_ZOOM <- 19
ESRI_GRAY_MAX_NATIVE_ZOOM    <- 16
MAP_MAX_ZOOM                 <- 19

# NOTE: do not set leafletOptions(preferCanvas = TRUE). L.Minichart extends
# L.CircleMarker, so on canvas the pie charts lose their DOM container and
# silently vanish (120 chart layers drop to 1).

# ================================================================================
# DATA — built by Tools/build_app_cache.R (runs once per R process)
# ================================================================================
DATA_DIR <- "data"

local({
  need   <- file.path(DATA_DIR, c("genes.rds", "map_layers.rds"))
  absent <- need[!file.exists(need)]
  if (length(absent)) {
    stop("Missing ", paste(absent, collapse = ", "), ". From the repository root run:\n",
         "  Rscript Tools/build_app_cache.R\n",
         "and start the app with Dashboard/ as the working directory.", call. = FALSE)
  }
})

# genes: one entry per folder in Input_Files/genes/ — map table, haplotype
# palette, years, StrainHub + country networks, tree figure path.
genes      <- readRDS(file.path(DATA_DIR, "genes.rds"))
map_layers <- readRDS(file.path(DATA_DIR, "map_layers.rds"))
shipping   <- map_layers$shipping
geno       <- map_layers$genomics   # NULL when there is no Genomics_stats.xlsx
ddrad      <- geno$ddrad_sites

# Same rule as the build: "Saudi Arabia" == "SaudiArabia".
norm_country <- function(x) tolower(gsub("[^A-Za-z]", "", x))

if (!seed_pubmed_cache(file.path(DATA_DIR, "pubmed_snapshot.rds"))) {
  message("No PubMed snapshot found — the Literature tab will fetch live on first visit.")
}

# Shipping routes are stored as one pre-serialised GeoJSON Feature per route.
# Selecting routes is then string pasting, not geometry serialisation.
ship_geojson <- function(ids) {
  paste0('{"type":"FeatureCollection","features":[',
         paste(shipping$features[ids], collapse = ","), "]}")
}

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

# leaflet.minicharts ships its JavaScript attached to the addMinicharts() call.
# Over a leafletProxy that races: the method is invoked before the script has
# loaded and is dropped as "Unknown method addMinicharts". Attaching the
# dependency to the initial widget loads it with the page instead.
minichart_deps <- utils::getFromNamespace("minichartDeps", "leaflet.minicharts")

create_base_map <- function(tile_size = 256) {
  m <- leaflet(options = leafletOptions(maxZoom = MAP_MAX_ZOOM)) %>%
    # Below the default overlay pane (400) that the pie charts draw in, so no
    # overlay line is ever drawn over a pie.
    addMapPane("countryStatus", zIndex = 200) %>%
    addMapPane("shipping",      zIndex = 250) %>%
    addMapPane("railways",      zIndex = 300) %>%
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
    # Railways start visible but are pushed in by the server once the map has
    # rendered; Shipping starts hidden and is only sent if ticked.
    hideGroup("Shipping Routes")

  m$dependencies <- c(m$dependencies, minichart_deps())
  m
}

# ---- Genomic layers (WGS / ddRAD circles) --------------------------------------
# WGS (per country) and ddRAD (per site) stay separate layers on purpose: ddRAD
# reads a small fraction of the genome, so its values aren't comparable to WGS.
geno_tip   <- function(key) geno$descriptors$tip[geno$descriptors$key == key]
fmt_num    <- function(x, digits) ifelse(is.na(x), "—", formatC(x, format = "f", digits = digits))
split_refs <- function(x) unique(trimws(unlist(strsplit(x, ";\\s*(?=[A-Z])", perl = TRUE))))

# One pop-up table row; the descriptor sheet's tooltip rides on `key`.
stat_row <- function(label, value, key = "") {
  tip <- if (key %in% geno$descriptors$key) geno_tip(key) else key   # a descriptor key, or the tooltip text itself
  sprintf("<tr title=\"%s\"><td style='padding:1px 10px 1px 0;color:#555'>%s</td><td>%s</td></tr>",
          htmltools::htmlEscape(tip, attribute = TRUE), label, value)
}

stat_table <- function(...) paste0("<table style='font-size:12px;border-collapse:collapse'>",
                                   paste0(c(...), collapse = ""), "</table>")

# Pop-up: place name, a subtitle, the stat tables, a source line with a link to
# the About page's explanation of the statistics.
geno_popup <- function(title, sub, body, src) paste0(
  "<div style='font-family:sans-serif;min-width:170px'><b style='font-size:13px'>",
  htmltools::htmlEscape(title), "</b>",
  "<div style='color:#666;font-size:11px;margin:3px 0 6px'>", sub, "</div>", body,
  "<div style='color:#888;font-size:11px;margin-top:5px'>", src,
  " &middot; <a href='#' onclick=\"Shiny.setInputValue('show_genomic_about', Date.now()); return false;\">",
  "What do these mean?</a></div></div>")

# One circle per place: data.frame(id, lat, lon, pi, popup).
#  ddrad: the sheet's own sites, each with its site values and its region's.
#  wgs:   WGS is per country, so one circle per country.
geno_points <- function(kind) {
  if (kind == "ddrad") {
    s   <- geno$ddrad_sites
    reg <- geno$ddrad_regions
    popup <- vapply(seq_len(nrow(s)), function(i) {
      r <- if (is.null(reg)) NULL else reg[match(s$region[i], reg$region), ]
      src <- htmltools::htmlEscape(s$short[i])
      if (!is.na(s$doi[i])) src <- sprintf("<a href='%s' target='_blank' rel='noopener'>%s</a>", s$doi[i], src)
      geno_popup(s$site[i], paste0("ddRAD-seq &middot; ", htmltools::htmlEscape(s$country[i])),
        paste0(
          stat_table(stat_row("Samples", s$n[i], "Individuals genotyped at this site"),
                     stat_row("&pi;", fmt_num(s$pi[i], 3), "pi"),
                     stat_row("S", format(s$S[i], big.mark = ","), "Segregating (polymorphic) sites"),
                     stat_row("Tajima&#39;s D", fmt_num(s$tajima_d[i], 2), "tajima_d")),
          if (!is.null(r) && !is.na(r$region)) paste0(
            "<div style='margin-top:6px;font-size:12px;font-weight:600'>", htmltools::htmlEscape(r$region),
            " <span style='font-weight:400;color:#666'>(", r$sites, " sites, ", r$n, " samples)</span></div>",
            stat_table(stat_row("&pi;", fmt_num(r$pi, 3), "pi"),
                       stat_row("Tajima&#39;s D", fmt_num(r$tajima_d, 2), "tajima_d"),
                       stat_row("Fst", fmt_num(pmax(r$fst, 0), 3), "fst")))
        ), src)
    }, "")
    return(data.frame(id = paste0("ddrad_", seq_len(nrow(s))), lat = s$lat, lon = s$lon, pi = s$pi, popup = popup))
  }

  # One circle per WGS country, on a point inside its polygon.
  w  <- geno$wgs
  cs <- map_layers$country_status
  i  <- match(norm_country(w$country), norm_country(cs$admin))
  xy <- sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(sf::st_geometry(cs)[i])))
  popup <- vapply(seq_len(nrow(w)), function(j) geno_popup(
    w$country[j], "Whole-genome sequencing (WGS) &middot; country-wide values",
    stat_table(stat_row("Samples", w$n[j], "Mosquitoes sequenced in this country"),
               stat_row("Sites", w$sites[j], "Collection sites sequenced in this country"),
               stat_row("&pi;", fmt_num(w$pi[j], 4), "pi"),
               stat_row("Tajima&#39;s D", fmt_num(w$tajima_d[j], 2), "tajima_d"),
               stat_row("Fst", fmt_num(pmax(w$fst[j], 0), 3), "fst")),
    htmltools::htmlEscape(paste(split_refs(w$source[j]), collapse = "; "))), "")
  data.frame(id = paste0("wgs_", seq_len(nrow(w))), lat = xy[, "Y"], lon = xy[, "X"], pi = w$pi, popup = popup)
}

# Circles coloured on a yellow-to-red heat scale of nucleotide diversity. Each
# layer gets its own scale: WGS and ddRAD values are ~10x apart.
add_pi_circles <- function(proxy, pts, digits) {
  if (is.null(pts)) return(proxy)
  pal <- colorNumeric("YlOrRd", domain = range(pts$pi, na.rm = TRUE))
  proxy %>%
    addCircleMarkers(
      lng = pts$lon, lat = pts$lat, layerId = pts$id, group = "Genomic circles",
      radius = 11, weight = 1.5, color = "#444", fillColor = pal(pts$pi), fillOpacity = 0.9,
      popup = pts$popup
    ) %>%
    addLegend(position = "bottomleft", pal = pal, values = pts$pi, opacity = 1,
              title = "Nucleotide<br>diversity (&pi;)", layerId = "pi_legend",
              labFormat = labelFormat(digits = digits))
}

# Description card text for the tree and network tabs.
explainer <- function(description, how, glossary) tagList(
  tags$p(class = "help-text", description),
  tags$h6(class = "explainer-head", "How to read it"),
  tags$div(class = "help-text", how),
  tags$details(
    class = "explainer-glossary",
    tags$summary("Glossary"),
    tags$dl(class = "help-text", lapply(names(glossary), function(k)
      tagList(tags$dt(k), tags$dd(glossary[[k]]))))
  )
)

# Descriptions of the three genomic statistics, worded for the WGS analysis
# (shown in the map sidebar and on the Genomic Statistics page). Hover text on
# the stats comes from the workbook's "Tool Tip" column instead.
GENO_TEXT <- list(
  pi = list(
    title = "Nucleotide diversity (π)",
    description = paste(
      "Nucleotide diversity. This measure shows how genetically varied the mosquitoes within each",
      "population are. The genome was divided into 10,000-letter stretches (windows). In each window,",
      "DNA was compared between every pair of mosquitoes from the same population and the proportion",
      "of differing letters was calculated. These values were then averaged across windows to give one",
      "number per population."),
    how = tags$p(paste(
      "Higher values mean mosquitoes from that country carry more varied DNA, as you would expect from",
      "a large or long-established population. Lower values mean the mosquitoes are more genetically",
      "alike, which can happen when a population descends from a small number of founders, such as after",
      "a recent introduction. Values are small numbers because most of the genome is identical between",
      "any two mosquitoes.")),
    glossary = c(
      "Nucleotide" = "One letter (A, T, C, or G) in the DNA code.",
      "Variant" = "A spot in the genome where mosquitoes differ from one another.",
      "Nucleotide diversity (π)" = "The average proportion of DNA letters that differ between two randomly chosen mosquitoes from the same group.",
      "Founder effect" = "Reduced variation in a new population started by a small number of individuals.")
  ),

  tajima_d = list(
    title = "Tajima's D",
    description = paste(
      "Tajima's D. This statistic compares two ways of measuring genetic variation to reveal clues about",
      "a population's recent history. It was calculated in 10,000-letter windows and averaged across the",
      "genome for each population. Populations with fewer samples give less stable values, so sample",
      "sizes are shown alongside."),
    how = tagList(
      tags$ul(
        tags$li(tags$b("Near zero:"), " the population looks stable, with a normal balance of common and rare variants."),
        tags$li(tags$b("Positive:"), " rare variants are scarce. This can follow a population crash (bottleneck) or the mixing of genetically distinct groups."),
        tags$li(tags$b("Negative:"), " rare variants are abundant. This can follow a recent expansion, when a population grows quickly and new variants appear.")
      ),
      tags$p("These are clues about history, not proof. Sampling choices and other forces can push the value up or down.")
    ),
    glossary = c(
      "Rare variant" = "A DNA difference found in only one or a few mosquitoes in a given population.",
      "Bottleneck" = "A sharp drop in population size that removes genetic variation.",
      "Population expansion" = "Rapid growth of a population, often after reaching new territory.",
      "Tajima's D" = "A test that compares the number of variable sites with the average pairwise difference. The gap between them signals whether the population recently shrank or grew.",
      "Window" = "A fixed-length stretch of the genome (here 10,000 or 100,000 letters) analyzed as a unit so patterns can be compared along the chromosome.")
  ),

  fst = list(
    title = "Fst",
    description = paste(
      "Fst. This value shows how genetically distinct two groups of mosquitoes are. Each country's",
      "mosquitoes were divided into two groups based on geography, except for Djibouti, where all",
      "collections came from Djibouti City and were therefore divided randomly into two groups. Fst was",
      "calculated across the whole genome, and plots show the same comparison in 100,000-letter windows."),
    how = tags$p(
      "Values run from 0 to 1. Values near 0 mean the two groups are essentially one interbreeding
      population, with mosquitoes moving freely between them. Higher values mean the groups have drifted
      apart, which suggests limited mixing, for example because of distance or barriers. Random splits
      are expected to give values near 0, so they act as a baseline for comparison. Negative values
      should be interpreted as 0. These values compare two groups of mosquitoes from within the same
      country, so they only measure local structure. They should not be read as a measure of how
      different ", tags$em("An. stephensi"), " populations are overall, since the comparison here
      evaluates within-country differences."),
    glossary = c(
      "Population structure" = "Genetic differences that divide a species into partly separate groups.",
      "Gene flow" = "The movement and mixing of genes between groups when individuals breed across them.",
      "Fst" = "The share of total genetic variation that is explained by differences between groups rather than within them.",
      "Weighted Fst" = "A genome-wide value that combines all variants by adding up their evidence, so informative sites count for more than uninformative ones.",
      "Window" = "A fixed-length stretch of the genome (here 10,000 or 100,000 letters) analyzed as a unit so patterns can be compared along the chromosome.")
  )
)

add_pie_markers <- function(map, data, color_palette,
                            lat_col       = "Latitude",
                            lon_col       = "Longitude",
                            category_col  = "Haplotype",
                            city_col      = "City.Town",
                            state_col     = "County.State",
                            country_col   = "Country",
                            radius        = 20) {

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
      country = {
        countries <- .data[[country_col]][!is.na(.data[[country_col]]) & .data[[country_col]] != ""]
        if (length(countries) > 0) names(which.max(table(countries))) else NA_character_
      },
      .groups = "drop"
    )

  df_pie <- data %>%
    filter(!is.na(.data[[lat_col]]), !is.na(.data[[lon_col]])) %>%
    count(.data[[lat_col]], .data[[lon_col]], .data[[category_col]]) %>%
    pivot_wider(names_from = all_of(category_col), values_from = n, values_fill = 0) %>%
    left_join(city_lookup, by = c(lat_col, lon_col))

  coords <- df_pie %>% select(all_of(c(lat_col, lon_col)))
  pies   <- df_pie %>% select(-all_of(c(lat_col, lon_col, "city", "country")))

  # The palette is built in haplotype-number order, so following it orders the
  # slices the same way and guarantees every slice has a colour.
  pies   <- pies[, intersect(names(color_palette), colnames(pies))]
  colors <- unname(color_palette[colnames(pies)])

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
      rows_html, "</table>",
      "</div>"
    )
  }

  popup_html <- unlist(popup_list[layer_ids])

  map %>%
    addMinicharts(
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
    ) %>%
    # Invisible markers over each pie: their clicks arrive as input$map_marker_click,
    # which drives the shipping-route highlight.
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
}

add_country_status_layer <- function(map, world = map_layers$country_status) {
  pal <- colorFactor(
    palette  = c("native" = "#4CAF50", "invasive" = "#FF9800"),
    domain   = c("native", "invasive"),
    na.color = "transparent"
  )

  map %>%
    addPolygons(
      data        = world,
      fillColor   = ~pal(status_clean),
      fillOpacity = 0.15,
      color       = "#444444",
      weight      = 0.5,
      opacity     = 0.5,
      group       = "Country Status",
      options     = pathOptions(pane = "countryStatus"),
      label       = ~paste0(admin, ": ", ifelse(is.na(status_clean), "No Data", status_clean))
    ) %>%
    addLegend(
      position = "bottomright",
      pal      = pal,
      values   = c("native", "invasive"),
      title    = "Native vs. Invasive",
      group    = "Country Status"
    )
}

# Pre-serialised at build time: the browser JSON.parse()s it, R does no geometry
# work. Drawing this from sf used to cost 2.6 s of R time per visitor.
add_railway_layer <- function(map) {
  types <- map_layers$railway_types
  map %>%
    addGeoJSON(
      map_layers$railways_geojson,
      group   = "Railways",
      weight  = 1.2,
      opacity = 0.75,
      fill    = FALSE,
      options = pathOptions(pane = "railways")
    ) %>%
    addLegend(
      position = "bottomright",
      colors   = unname(map_layers$railway_colours[types]),
      labels   = types,
      title    = "Railway type",
      group    = "Railways",
      layerId  = "railway_legend"
    )
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
  success       = "#4CAF50",
  warning       = "#FF7F00",
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
  # Clicking the logo or name goes back to the map.
  title = tags$a(
    href    = "#",
    class   = "brand-home",
    title   = "Back to the map",
    onclick = "document.querySelector('.navbar a[data-value=\"map_tab\"]').click(); return false;",
    tags$img(src = "logo.png", height = "100px", alt = "", style = "margin-right:10px;vertical-align:middle"),
    "The Vector Invasion Atlas"
  ),
  theme    = app_theme,
  id       = "main_nav",
  fillable = FALSE,

  header = tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
    tags$link(rel = "icon", type = "image/png", href = "favicon.png"),
    tags$script(HTML(
      "Shiny.addCustomMessageHandler('scroll_to', function(id) {
         setTimeout(function() { var el = document.getElementById(id);
           if (el) el.scrollIntoView({behavior: 'smooth'}); }, 150);
       });

       // Pan/zoom for the tree figure. The SVG is inlined and zoomed by editing
       // its viewBox, so it stays sharp at any zoom and only the visible part is drawn.
       function treeZoom(boxId, src) {
         var box = document.getElementById(boxId);
         fetch(src).then(function(r) { return r.text(); }).then(function(txt) {
           var view = box.querySelector('.tree-zoom-view');
           view.innerHTML = txt;
           var svg = view.querySelector('svg');
           var w = +svg.getAttribute('width'), h = +svg.getAttribute('height');
           if (!svg.getAttribute('viewBox')) svg.setAttribute('viewBox', '0 0 ' + w + ' ' + h);
           var home = svg.getAttribute('viewBox').split(/[ ,]+/).map(Number), vb = home.slice();
           svg.removeAttribute('width'); svg.removeAttribute('height');
           function set(v) { vb = v; svg.setAttribute('viewBox', v.join(' ')); }
           function toSvg(cx, cy) {                       // screen -> SVG coordinates
             var p = svg.createSVGPoint(); p.x = cx; p.y = cy;
             return p.matrixTransform(svg.getScreenCTM().inverse());
           }
           function zoom(f, cx, cy) {
             f = Math.min(f, vb[2] / (home[2] / 40)); f = Math.max(f, vb[2] / home[2]);   // 1x-40x
             var q = toSvg(cx, cy);
             set([q.x - (q.x - vb[0]) / f, q.y - (q.y - vb[1]) / f, vb[2] / f, vb[3] / f]);
           }
           function centre() { var r = view.getBoundingClientRect(); return [r.left + r.width / 2, r.top + r.height / 2]; }
           view.addEventListener('wheel', function(e) {
             e.preventDefault(); zoom(e.deltaY < 0 ? 1.25 : 0.8, e.clientX, e.clientY);
           }, {passive: false});
           var last = null;
           view.addEventListener('pointerdown', function(e) {
             last = [e.clientX, e.clientY]; view.setPointerCapture(e.pointerId); view.classList.add('dragging');
           });
           view.addEventListener('pointermove', function(e) {
             if (!last) return;
             var k = svg.getScreenCTM().a;               // screen px per SVG unit
             set([vb[0] - (e.clientX - last[0]) / k, vb[1] - (e.clientY - last[1]) / k, vb[2], vb[3]]);
             last = [e.clientX, e.clientY];
           });
           function stop() { last = null; view.classList.remove('dragging'); }
           view.addEventListener('pointerup', stop); view.addEventListener('pointercancel', stop);
           box.querySelectorAll('[data-zoom]').forEach(function(b) {
             b.addEventListener('click', function() {
               var c = centre(), z = b.getAttribute('data-zoom');
               if (z === 'reset') set(home.slice()); else zoom(+z, c[0], c[1]);
             });
           });
         });
       }"))
  ),

  # ---- Gene selector (top-right of navbar) -----------------------------------
  # One entry per folder in Input_Files/genes/, after Tools/build_app_cache.R.
  nav_item(
    tags$div(
      style = "padding: 6px 10px;",
      selectInput(
        "selected_gene",
        label    = NULL,
        choices  = names(genes),
        selected = names(genes)[1],
        width    = "120px"
      )
    )
  ),

  # ---- Tab 1: Interactive Map ------------------------------------------------
  nav_panel(
    title = tags$span(icon("globe"), "Distribution Map"),
    value = "map_tab",

    tags$div(
      class = "map-banner",
      tags$h2("Mapping the invasion of ", tags$em("Anopheles stephensi")),
      tags$p(
        "An urban malaria vector native to South Asia and the Arabian Peninsula, ",
        tags$em("An. stephensi"), " has spread across the Horn of Africa and beyond ",
        "since 2012. This atlas brings together COI mitochondrial haplotypes from ",
        "field collections to show where each lineage has been found, and which ",
        "populations are genetically connected."
      ),
      tags$div(
        class = "banner-hint",
        icon("circle-info"),
        tags$span(
          " Each pie chart is one collection site, split by haplotype — click it for a breakdown. ",
          "Use the ", tags$strong("Temporal Filter"), " to watch the range expand through time, ",
          "and the ", tags$strong("layers control"), " on the map to overlay railways and ",
          "shipping routes as possible corridors of spread."
        )
      )
    ),

    layout_sidebar(
      fillable = FALSE,

      sidebar = sidebar(
        width = 280,
        open  = TRUE,

        conditionalPanel(
          condition = "!input.map_data || input.map_data === 'haplotype'",
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

        if (!is.null(geno)) conditionalPanel(
          condition = "input.map_data === 'wgs' || input.map_data === 'ddrad'",
          tags$div(class = "sidebar-section-header", "Genomic statistics"),
          uiOutput("genomic_info"),
          tags$details(
            class = "geno-key",
            tags$summary("How to read these statistics"),
            lapply(GENO_TEXT, function(g) tags$details(
              class = "geno-stat",
              tags$summary(g$title),
              explainer(g$description, g$how, g$glossary)
            )),
            tags$p(tags$a(href = "#", onclick = "Shiny.setInputValue('show_genomic_about', Date.now()); return false;",
                          "Sources and more: About › Genomic Statistics"))
          )
        )
      ),

      card(
        full_screen = TRUE,
        class       = "map-card",
        if (!is.null(geno)) card_header(
          class = "map-data-header",
          tags$label(`for` = "map_data", "Show on map"),
          selectInput("map_data", NULL, selectize = FALSE, width = "130px",
                      choices = c("Haplotype" = "haplotype", "WGS" = "wgs", "ddRAD" = "ddrad"))
        ),
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

      card(
        card_header(tags$span(icon("circle-info"), " About this Tree")),
        class = "controls-card",
        explainer(
          tagList(
            "Maximum likelihood phylogenetic tree of ", tags$em("Anopheles stephensi"), " COI
            sequences. This family tree shows how mosquitoes from different places are related,
            based on one short stretch of mitochondrial DNA (the COI gene). It was inferred with
            RAxML and drawn with iTOL. Scroll or use the + and − buttons to zoom, drag to move
            around, and use Download SVG to save the full figure."
          ),
          tags$p(
            "Each tip (end of a branch) is one sequenced mosquito, labelled with its sample ID and
            country. Up to 20 mosquitoes per country and year are included, so heavily sampled
            places don't dominate the tree. Tips that sit close together have more similar
            sequences, and mosquitoes that share a haplotype sit side by side. Groups of tips that
            share a common branch point are called clades; think of them as extended families of
            related sequences. Longer branches mean more DNA differences between sequences.
            Numbers at branch points show how strongly the data support that grouping, with higher
            values meaning more confidence."
          ),
          list(
            "Phylogenetic tree" = "A diagram of evolutionary relationships, like a family tree.
              Sequences that share a more recent common ancestor sit closer together.",
            "COI" = "A gene in mosquito mitochondria (the cell's energy-producing structures) that
              is read in many species. It changes slowly enough to compare populations and is
              inherited only from the mother, so the tree traces maternal lineages.",
            "Haplotype" = "One distinct version of a DNA sequence. Many mosquitoes can share the
              same haplotype.",
            "Clade" = "A branch and everything that descends from it: a group of sequences that
              likely share a common ancestor.",
            "Branch length" = "The amount of genetic change between points on the tree. It
              reflects DNA differences, not distance in kilometres or years.",
            "Maximum likelihood" = "A statistical method that tests many possible trees and picks
              the one that best explains the observed DNA differences.",
            "Bootstrap support" = "The numbers at branch points: a confidence score from repeating
              the analysis on resampled data. Higher values mean a grouping appeared more
              consistently.",
            "iTOL" = "Interactive Tree Of Life, the web tool used to draw this figure."
          )
        )
      ),

      card(
        full_screen = TRUE,
        uiOutput("tree_view")
      )
    )
  ),

  # ---- Tab 3: Gene Network ----------------------------------------------------
  nav_panel(
    title = tags$span(icon("circle-nodes"), "Gene Network"),
    value = "network_tab",

    layout_columns(
      col_widths = c(3, 9),

      # Wording follows the reviewer responses: arrows are inferred country changes
      # on the COI tree, not observed movement, and the betweenness ranking is not
      # stable under resampling, so no country is called a hub.
      card(
        card_header(tags$span(icon("circle-info"), " About this Network")),
        class = "controls-card",
        explainer(
          tagList(
            "This network uses one mitochondrial gene (COI) to estimate which countries ",
            tags$em("An. stephensi"), " lineages may have moved between. Each circle is a country.
            Each arrow points from the country inferred for a lineage's ancestor on the COI tree to
            the country where its descendants were sampled, and thicker arrows mean that change was
            inferred more often. Larger circles have higher betweenness centrality: they lie
            between many other countries in this network. These are hypotheses drawn from genetic
            patterns, not observed movements of mosquitoes."
          ),
          tags$ul(
            tags$li(tags$b("Circles: "), "many of the network's shortest routes between other
                    countries pass through this country."),
            tags$li(tags$b("Arrows: "), "the tree repeatedly places an ancestor in one country
                    and its descendants in the other, in that direction. Hover over an arrow for
                    the count."),
            tags$li(tags$b("Arrow direction: "), "from the inferred ancestral country toward the
                    country where descendants were sampled."),
            tags$li(tags$b("No arrow: "), "no change between those countries was inferred in this
                    dataset. It does not prove there was none: countries with few samples have
                    fewer chances to show links.")
          ),
          list(
            "Node" = "A country in the dataset.",
            "Arrow" = "An inferred change of country between an ancestor and its descendants on the
              COI tree: a possible dispersal pathway, not a measured movement of mosquitoes.",
            "Arrow thickness" = "How many times that change was inferred.",
            "Betweenness centrality" = "How often a country lies on the shortest path between two
              other countries in this network. It reflects the network's shape, which depends on
              which countries were sampled and how heavily. It does not measure how many mosquitoes
              live there.",
            "Ancestral state" = "The country inferred for a common ancestor on the family tree,
              used to decide where a lineage came from.",
            "Maximum parsimony" = "A method that assigns each ancestor the country that explains the
              data with the fewest moves between countries.",
            "COI" = "A mitochondrial gene inherited from the mother, so this network traces maternal
              lineages only."
          )
        )
      ),

      card(
        full_screen = TRUE,
        card_header("StrainHub: Phylogeographic Dispersal Network"),
        withSpinner(
          visNetworkOutput("haplo_network", height = "700px"),
          color = "#1a4f72"
        )
      )
    ),

    layout_columns(
      col_widths = c(3, 9),

      card(
        hr(class = "sidebar-divider"),
        tags$p(class = "help-text",
               style = "padding: 72px; font-size: 16px; margin: 8px 0;",
               icon("circle-info"),
               " Node size reflects number of included samples (capped at 20 per country/year). Edge thickness reflects the number of shared unique haplotypes between countries.",
               " Countries with more samples have more chances to detect shared haplotypes, so thin sampling can hide links.")
      ),

      card(
        full_screen = TRUE,
        card_header("An. stephensi: Country Haplotype-Sharing Network"),
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
      value = "about_tab",
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
          tags$p("For questions or data requests, contact: ",
                 tags$a(href = "mailto:Tamar_Carter@baylor.edu", "Tamar_Carter@baylor.edu")),
          tags$h5("Citation"),
          tags$p(tags$em("Coming Soon")),
          tags$hr(),
          tags$p(class = "text-muted small",
                 paste0("Built with R Shiny & bslib | Last updated: ", format(Sys.Date(), "%B %Y")))
        )
      )
    ),

    # Opened by the "What do these mean?" link in a genomic pop-up and the one in the sidebar.
    nav_panel(
      "Genomic Statistics",
      value = "genomics_tab",
      card(
        class = "about-card",
        card_body(
          tags$h4(id = "genomics-top", "Genomic statistics (WGS and ddRAD)"),
          if (is.null(geno)) tags$p(
            "No genomic statistics yet. Add ", tags$code("Input_Files/Genomics_stats.xlsx"),
            " and rerun ", tags$code("Rscript Tools/build_app_cache.R"), "."
          ) else tagList(
            tags$p(paste0(
              "Use the ", "Show on map", " menu above the map to switch from haplotype pie ",
              "charts to genome-wide nucleotide diversity (π), drawn as circles on a yellow-to-red ",
              "scale. Click a circle for the other statistics. The two kinds of sequencing are kept ",
              "apart because their values can't be compared. Whole-genome sequencing (WGS) is ",
              "summarized per country, so it appears as one circle per country. Double-digest restriction-site-associated DNA sequencing (ddRAD) is ",
              "reported per site, with the site's region alongside, for ", nrow(ddrad), " sites in ",
              paste(sort(unique(ddrad$country)), collapse = ", "), ". ", geno$ddrad_note, " ",
              "These come from the published studies below, not from the COI data, and don't ",
              "change with the year filter."
            )),
            tags$h5("How to read these statistics"),
            lapply(GENO_TEXT, function(g) tagList(tags$h6(class = "geno-stat-title", g$title),
                                                  explainer(g$description, g$how, g$glossary))),
            tags$p(
              "These descriptions are written for the WGS analysis. ddRAD reports the same ",
              "statistics from a much smaller part of the genome, so its values aren't comparable ",
              "to WGS. Two more numbers appear in the ddRAD pop-ups:"
            ),
            tags$ul(
              tags$li(tags$b("S (segregating sites): "),
                      "the number of variable DNA positions found at the site. S rises with the
                      number of mosquitoes sequenced, so compare sites with similar sample sizes."),
              tags$li(tags$b("Samples: "), "the number of mosquitoes sequenced.")
            ),
            tags$h5("Sources"),
            tags$ul(
              lapply(split_refs(geno$wgs$source), function(s) tags$li(s, " (WGS)")),
              lapply(split_refs(c(ddrad$source, geno$ddrad_regions$source)), function(s) {
                doi <- regmatches(s, regexpr("https?://doi\\.org/\\S+[^.\\s]", s, perl = TRUE))
                if (!length(doi)) return(tags$li(s, " (ddRAD)"))
                tags$li(sub(doi, "", s, fixed = TRUE),
                        tags$a(href = doi, target = "_blank", rel = "noopener", doi), " (ddRAD)")
              })
            )
          )
        )
      )
    )
  )
)

# ================================================================================
# SERVER
# ================================================================================

server <- function(input, output, session) {

  gene <- reactive(genes[[req(input$selected_gene)]])

  # Year sliders follow the selected gene's data.
  observeEvent(gene(), {
    yrs <- gene()$years
    updateSliderInput(session, "year_filter", min = min(yrs), max = max(yrs), value = max(yrs))
    updateSliderInput(session, "year_range",  min = min(yrs), max = max(yrs), value = range(yrs))
  })

  # Track selected location for shipping route highlighting
  selected_location <- reactiveVal(NULL)

  # ==============================================================================
  # TAB 1 — MAP
  # ==============================================================================

  filtered_data <- reactive({
    d <- gene()$haplo
    switch(input$time_mode,
           "exact"      = d %>% filter(Year == input$year_filter),
           "cumulative" = d %>% filter(Year <= input$year_filter),
           "range"      = d %>% filter(Year >= input$year_range[1],
                                       Year <= input$year_range[2]),
           "all"        = d
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

  # The map itself renders once per session; everything gene- or filter-dependent
  # is pushed through leafletProxy so the map keeps its zoom and layers.
  output$map <- renderLeaflet({
    create_base_map() %>% add_country_status_layer()
  })

  # "Show on map": haplotype pies (default), or WGS / ddRAD circles coloured by pi.
  map_kind <- reactive(if (is.null(input$map_data)) "haplotype" else input$map_data)

  observe({
    proxy <- leafletProxy("map") %>% removeControl("haplotype_legend")
    if (map_kind() != "haplotype") return()
    pal <- gene()$palette
    proxy %>%
      addLegend(
        position = "bottomleft",
        colors   = unname(pal),
        labels   = names(pal),
        title    = gene()$marker_label,   # "Haplotypes", "P47Rec Genotype", ...
        opacity  = 1,
        layerId  = "haplotype_legend"
      )
  })

  # ---- Overlays ---------------------------------------------------------------
  # input$map_groups holds the overlay groups currently ticked. Geometry is never
  # cleared on untick — Leaflet hides the group client-side — so switching a
  # layer back on is instant and re-sends nothing.
  railways_drawn     <- reactiveVal(FALSE)
  shipping_drawn_ids <- reactiveVal(NULL)

  layer_on <- function(name) isTRUE(name %in% input$map_groups)

  observeEvent(input$map_groups, {
    if (layer_on("Railways") && !railways_drawn()) {
      leafletProxy("map") %>% add_railway_layer()
      railways_drawn(TRUE)   # static geometry: at most once per session
    }
  }, ignoreNULL = TRUE)

  # Pie charts follow the gene and the temporal filter.
  observe({
    kind  <- map_kind()
    fd    <- filtered_data()
    proxy <- leafletProxy("map") %>% clearMinicharts() %>% clearMarkers()
    if (kind == "haplotype" && nrow(fd) > 0) {
      add_pie_markers(proxy, fd, color_palette = gene()$palette, radius = 30)
    }
  })

  # WGS / ddRAD circles follow neither the gene nor the temporal filter: the
  # genomic statistics come from their own workbook and aren't split by year.
  observe({
    kind  <- map_kind()
    proxy <- leafletProxy("map") %>% clearGroup("Genomic circles") %>% removeControl("pi_legend")
    if (kind != "haplotype" && !is.null(geno)) {
      add_pi_circles(proxy, geno_points(kind), digits = if (kind == "wgs") 4 else 3)
    }
  })

  # Handle pie chart clicks via invisible markers
  observeEvent(input$map_marker_click, {
    id <- input$map_marker_click$id
    if (is.null(id)) return(selected_location(NULL))
    country <- gene()$city_to_country[id]
    if (!is.na(country)) selected_location(unname(country))
  })

  # Sidebar blurb for the selected genomic layer.
  output$genomic_info <- renderUI({
    switch(map_kind(),
      wgs = tags$p(class = "geno-hint",
        "Whole-genome data, summarized per country: one circle per sequenced country. ",
        "Click a circle for all of its statistics."),
      ddrad = tags$p(class = "geno-hint",
        "Reduced-representation data, one circle per ddRAD site. Click a circle for the ",
        "site's statistics and its region's. ", geno$ddrad_note),
      NULL)
  })

  # The link in a genomic pop-up opens the Genomic Statistics page.
  observeEvent(input$show_genomic_about, {
    nav_select("main_nav", "genomics_tab")
    session$sendCustomMessage("scroll_to", "genomics-top")
  })

  # Reset on empty map click
  observeEvent(input$map_click, {
    if (is.null(input$map_click$id)) selected_location(NULL)
  }, ignoreInit = TRUE)

  # Routes active under the current temporal filter.
  active_route_ids <- reactive({
    yrs <- shipping$years
    target <- switch(input$time_mode,
                     "exact"      = yrs[yrs == input$year_filter],
                     "cumulative" = yrs[yrs <= input$year_filter],
                     "range"      = yrs[yrs >= input$year_range[1] & yrs <= input$year_range[2]],
                     "all"        = yrs)
    sort(unique(shipping$route_years$route_id[shipping$route_years$year %in% target]))
  })

  # Redraw only while the layer is on, and only when the selected routes actually
  # change: in cumulative mode 17 of the 19 year steps select an identical set.
  observe({
    on  <- layer_on("Shipping Routes")
    ids <- active_route_ids()   # dependency taken even while hidden
    if (!on || identical(ids, shipping_drawn_ids())) return()
    shipping_drawn_ids(ids)

    proxy <- leafletProxy("map") %>% clearGroup("Shipping Routes")
    if (length(ids)) {
      proxy %>% addGeoJSON(ship_geojson(ids), group = "Shipping Routes",
                           color = "#0077b6", weight = 1, opacity = 0.4, fill = FALSE,
                           options = pathOptions(pane = "shipping"))
    }
  })

  # Highlight routes touching the selected country; added after the base routes
  # so it draws on top of them.
  observe({
    proxy <- leafletProxy("map") %>% clearGroup("Shipping Highlight")
    loc <- selected_location()
    if (!layer_on("Shipping Routes") || is.null(loc)) return()
    r   <- shipping$routes
    ids <- intersect(active_route_ids(), r$route_id[r$from == loc | r$to == loc])
    if (length(ids)) {
      proxy %>% addGeoJSON(ship_geojson(ids), group = "Shipping Highlight",
                           color = "#d62728", weight = 2.5, opacity = 0.85, fill = FALSE,
                           options = pathOptions(pane = "shipping"))
    }
  })

  # ==============================================================================
  # TAB 2 — PHYLOGENETIC TREE
  # ==============================================================================

  output$tree_view <- renderUI({
    name <- input$selected_gene
    svg  <- gene()$tree_svg
    if (is.null(svg)) {
      return(tagList(
        card_header(paste(name, "Maximum Likelihood Tree")),
        tags$p(class = "help-text", style = "padding: 24px;",
               "No tree figure for ", name, " yet. Put its iTOL export at ",
               tags$code(paste0("Input_Files/genes/", name, "/tree.svg")),
               " and rerun ", tags$code("Rscript Tools/build_app_cache.R"), ".")
      ))
    }
    tagList(
      card_header(
        paste(name, "Maximum Likelihood Tree"),
        tags$a(
          href     = svg,
          download = paste0("anopheles_stephensi_", name, "_tree.svg"),
          class    = "btn btn-sm btn-outline-primary float-end",
          icon("download"), " Download SVG"
        )
      ),
      tags$div(
        id = "tree_zoom", class = "tree-zoom",
        tags$div(class = "tree-zoom-view"),
        tags$div(
          class = "tree-zoom-btns",
          tags$button(type = "button", `data-zoom` = "1.5", title = "Zoom in", "+"),
          tags$button(type = "button", `data-zoom` = "0.667", title = "Zoom out", "−"),
          tags$button(type = "button", `data-zoom` = "reset", title = "Fit the whole tree", "Reset")
        )
      ),
      tags$script(HTML(sprintf("treeZoom('tree_zoom', '%s');", svg)))
    )
  })

  # ==============================================================================
  # TAB 3 — NETWORKS (both precomputed by Tools/build_app_cache.R)
  # ==============================================================================

  output$haplo_network <- renderVisNetwork({
    net <- gene()$strainhub
    validate(need(net, paste0("No StrainHub network for ", input$selected_gene,
                              ": add tree.newick to its folder in Input_Files/genes/ and rebuild.")))

    visNetwork(net$nodes, net$edges,
               main = "Phylogeographic dispersal network (COI)") %>%
      visNodes(shape = "ellipse") %>%
      visLayout(randomSeed = 12) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = list(enabled = TRUE)) %>%
      visPhysics(solver = "repulsion", repulsion = list(nodeDistance = 200))
  })

  output$country_network <- renderVisNetwork({
    net <- gene()$country_net
    w <- visNetwork(net$nodes, net$links, main = "Country haplotype-sharing network (COI)")

    # Node x/y were computed by visIgraphLayout() at build time. These are the
    # remaining settings visIgraphLayout() applies: the igraphlayout flag makes
    # vis.js rescale the normalised coordinates, and node physics stays off so
    # the layout holds. Keeps igraph off the server.
    w$x$igraphlayout <- list(type = "square")

    w %>%
      visNodes(physics = FALSE) %>%
      visLayout(randomSeed = 1) %>%
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
      visPhysics(solver = "repulsion", repulsion = list(nodeDistance = 300),
                 stabilization = FALSE)
  })

  # ==============================================================================
  # TAB 4 — LITERATURE (PubMed module)
  # ==============================================================================

  pubBrowserServer("pubs")
}

shinyApp(ui = ui, server = server)
