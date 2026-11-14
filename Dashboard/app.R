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
ddrad      <- map_layers$ddrad   # NULL when there is no ddRAD spreadsheet

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

# Genome-wide (ddRAD) block for the pop-up of a pie at lat/lon in `country`, or
# "" if no ddRAD site matches. A match is the nearest ddRAD site in the same
# country within DDRAD_MATCH_KM: the spreadsheet's coordinates for a town differ
# from the COI ones by up to ~3 km (Hargeisa), and Loyada (Djibouti) is 4 km
# from Lawyacado (Somalia), hence the country check.
DDRAD_MATCH_KM <- 5

ddrad_popup <- function(lat, lon, country) {
  if (is.null(ddrad) || !nrow(ddrad) || is.na(country)) return("")
  r  <- pi / 180
  km <- 6371 * 2 * asin(sqrt(sin((ddrad$lat - lat) * r / 2)^2 +
                             cos(lat * r) * cos(ddrad$lat * r) * sin((ddrad$lon - lon) * r / 2)^2))
  km[tolower(trimws(ddrad$country)) != tolower(trimws(country))] <- Inf
  i <- which.min(km)
  if (km[i] > DDRAD_MATCH_KM) return("")
  s <- ddrad[i, ]

  src <- htmltools::htmlEscape(s$short)
  if (!is.na(s$doi)) src <- sprintf("<a href='%s' target='_blank' rel='noopener'>%s</a>", s$doi, src)
  row <- function(label, value, tip) sprintf(
    "<tr title='%s'><td style='padding:1px 8px 1px 0;color:#555'>%s</td><td>%s</td></tr>", tip, label, value)

  paste0(
    "<div style='margin-top:8px;padding-top:6px;border-top:1px solid #e3e3e3'>",
    "<div style='font-size:12px;font-weight:600'>Genome-wide (ddRAD-seq) ",
    "<a href='#' title='What do these mean?' style='text-decoration:none' ",
    "onclick=\"Shiny.setInputValue('show_genomic_about', Date.now()); return false;\">&#9432;</a></div>",
    "<div style='color:#666;font-size:11px;margin:1px 0 4px'>Site: ", htmltools::htmlEscape(s$site), "</div>",
    "<table style='font-size:12px;border-collapse:collapse'>",
    row("Samples", s$n, "Individuals genotyped at this site"),
    row("&pi;", formatC(s$pi, format = "f", digits = 3), "Nucleotide diversity"),
    row("S", format(s$S, big.mark = ","), "Segregating (polymorphic) sites"),
    row("Tajima&#39;s D", formatC(s$tajima_d, format = "f", digits = 2), "Compares nucleotide diversity with the number of segregating sites"),
    "</table>",
    "<div style='color:#888;font-size:11px;margin-top:3px'>", src, "</div></div>"
  )
}

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
      ddrad_popup(coords[[lat_col]][i], coords[[lon_col]][i], df_pie$country[i]),
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
  title = tags$span(
    tags$img(src = "logo.png", height = "100px", style = "margin-right:10px;vertical-align:middle"),
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
       });"))
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

      card(
        card_header(tags$span(icon("circle-info"), " About this Tree")),
        class = "controls-card",
        tags$p(class = "help-text",
               "Maximum likelihood phylogenetic tree of", tags$em("Anopheles stephensi"),
               "haplotypes inferred with RAxML. Tree is hosted and visualised
              interactively via iTOL. Use the iTOL toolbar to zoom, collapse clades,
              and export the figure.")
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

      card(
        hr(class = "sidebar-divider"),
        tags$p(class = "help-text",
               style = "padding: 72px; font-size: 16px; margin: 8px 0;",
               icon("circle-info"),
               " Node size shows betweenness centrality: how often a country lies on the shortest path between other countries.",
               " Each arrow means a lineage sampled in one country has its closest inferred ancestor in another, based on the COI tree.",
               " Edge width counts how many such changes were inferred.",
               " Arrows are possible dispersal pathways from a single gene, not measured mosquito movement.")
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
          if (!is.null(ddrad) && nrow(ddrad)) tagList(
            tags$h5(id = "genomic-data", "Genome-wide data (ddRAD-seq)"),
            tags$p(paste0(
              "At ", nrow(ddrad), " collection sites in ",
              paste(sort(unique(ddrad$country)), collapse = ", "),
              ", the map pop-up also shows genome-wide statistics from double-digest ",
              "restriction-site-associated DNA sequencing (ddRAD-seq). They come from the ",
              "published studies below, not from the COI data, and don't change with the ",
              "year filter."
            )),
            tags$ul(
              tags$li(tags$b("Samples:"), " individuals genotyped at the site."),
              tags$li(tags$b("π (nucleotide diversity):"),
                      " the average genetic difference between two individuals from the site.
                      Higher values mean more genetic diversity."),
              tags$li(tags$b("S (segregating sites):"),
                      " the number of polymorphic positions found. S rises with sample size,
                      so compare sites with similar numbers of samples."),
              tags$li(tags$b("Tajima's D:"),
                      " compares π with S. Values near 0 fit a population of stable size;
                      positive values can reflect a recent bottleneck or population structure,
                      and negative values a recent expansion.")
            ),
            tags$p("Sources:"),
            tags$ul(lapply(unique(ddrad$source), function(s) {
              doi <- ddrad$doi[match(s, ddrad$source)]
              if (is.na(doi)) return(tags$li(s))
              tags$li(sub(doi, "", s, fixed = TRUE),
                      tags$a(href = doi, target = "_blank", rel = "noopener", doi))
            }))
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

  observeEvent(gene(), {
    pal <- gene()$palette
    leafletProxy("map") %>%
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
    fd    <- filtered_data()
    proxy <- leafletProxy("map") %>% clearMinicharts() %>% clearMarkers()
    if (nrow(fd) > 0) {
      add_pie_markers(proxy, fd, color_palette = gene()$palette, radius = 30)
    }
  })

  # Handle pie chart clicks via invisible markers
  observeEvent(input$map_marker_click, {
    id <- input$map_marker_click$id
    if (is.null(id)) return(selected_location(NULL))
    country <- gene()$city_to_country[id]
    if (!is.na(country)) selected_location(unname(country))
  })

  # The (i) in a pop-up's ddRAD block opens the About page's explanation.
  observeEvent(input$show_genomic_about, {
    nav_select("main_nav", "about_tab")
    session$sendCustomMessage("scroll_to", "genomic-data")
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
        paste0("iTOL — ", name, " Maximum Likelihood Tree"),
        tags$a(
          href     = svg,
          download = paste0("anopheles_stephensi_", name, "_tree.svg"),
          class    = "btn btn-sm btn-outline-primary float-end",
          icon("download"), " Download SVG"
        )
      ),
      tags$img(src = svg, width = "100%", style = "height:700px; object-fit:contain;")
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
