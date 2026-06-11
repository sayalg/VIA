# ================================================================================
# Haplotype Distribution Explorer — Multi-Tab Shiny App
# ================================================================================

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
library(readxl)
library(visNetwork)
library(DT)
library(ape)
library(ggplot2)
library(ggtree)
library(strainhub)
library(stringr)
library(treeio)
library(httr)
library(xml2)
library(plotly)

source("pubmed_module.R")

# ---- Configuration ----
MAPBOX_TOKEN      <- "pk.eyJ1Ijoic2F5YWxnIiwiYSI6ImNtaWdha3Y4eDA1YmczZXEybjZvZjE0YTQifQ.T_HqxBxbmcdViI2L0LIzgQ"
MAPBOX_STYLE      <- "https://api.mapbox.com/styles/v1/sayalg/cmkr2q0c3000q01s87umrehau/tiles/256/{z}/{x}/{y}@2x?access_token={accessToken}"
INVALID_LOCATIONS <- c("Laboratory")

# ================================================================================
# DATA LOADING
# ================================================================================

samples_df <- read.csv("../Input_Files/FinalCOI_metadata.csv", sep = ",", header = TRUE)
samples_df[samples_df == ""] <- NA

invasive_status <- read_excel("../Input_Files/MTM_INVASIVE_VECTOR_SPECIES_20251205.xlsx", sheet = "Data")
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

# Load pre-geocoded file (run fill_missing_coords() offline to regenerate)
# df_geo <- fill_missing_coords(samples_df)
# write.csv(df_geo, "../Input_Files_Latest/FinalCOI_metadata_w_coor.csv", row.names = FALSE)
df_geo   <- read.csv("../Input_Files/FinalCOI_metadata_w_coor.csv", sep = ",", header = TRUE)
df_haplo <- df_geo %>% filter(!is.na(Haplotype))

# ---- Cached world polygons (loaded once) ----
world_polygons <- ne_countries(returnclass = "sf") %>%
  mutate(name_clean = tolower(admin))

# ---- Africa railways (pre-built by prepare_africa_railways.R) ----
# Run prepare_africa_railways.R once offline to generate this file.
africa_railways_path <- "../Input_Files/africa_railways.rds"
africa_railways <- if (file.exists(africa_railways_path)) {
  readRDS(africa_railways_path)
} else {
  warning("africa_railways.rds not found — railway layer will be hidden. Run prepare_africa_railways.R.")
  NULL
}

# ---- Shipping routes (pre-built by prepare_shipping_routes.R) ----
# Named list: list("2006" = sf, "2007" = sf, ...)
# Each sf has columns: from, to, geometry (LINESTRING, WGS84)
shipping_rds_path <- "../Input_Files/shipping_routes.rds"
shipping_routes_by_year <- if (file.exists(shipping_rds_path)) {
  readRDS(shipping_rds_path)
} else {
  warning("shipping_routes.rds not found — shipping layer disabled. Run prepare_shipping_routes.R.")
  NULL
}
shipping_years <- if (!is.null(shipping_routes_by_year)) {
  as.integer(names(shipping_routes_by_year))
} else {
  integer(0)
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

# ================================================================================
# MAP FUNCTIONS
# ================================================================================

create_base_map <- function(tile_size = 256) {
  leaflet() %>%
    addMapPane("countryStatus", zIndex = 200) %>%
    addMapPane("charts",        zIndex = 400) %>%
    addProviderTiles(providers$CartoDB.Positron,  group = "Base") %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
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
    hideGroup("Shipping Routes")
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
      label       = ~paste0(
        "<b>", ifelse(is.na(name) | name == "", "Unnamed railway", name), "</b><br>",
        "Type: ", fclass
      ) %>% lapply(htmltools::HTML),
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

add_shipping_layer <- function(map, routes_sf) {
  map <- clearGroup(map, "Shipping Routes")
  if (is.null(routes_sf) || nrow(routes_sf) == 0) return(map)

  map %>%
    addMapPane("shipping", zIndex = 250) %>%
    addPolylines(
      data    = routes_sf,
      color   = "#0077b6",
      weight  = 1,
      opacity = 0.4,
      group   = "Shipping Routes",
      options = pathOptions(pane = "shipping")
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
    tags$img(src = "logo.png", height = "32px", style = "margin-right:10px;vertical-align:middle"),
    "AnStep Atlas"
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
      # Left: controls
      card(
        card_header(tags$span(icon("sliders"), " Tree Options")),
        class = "controls-card",
        selectInput(
          "tree_layout", "Layout:",
          choices  = c("Rectangular" = "rectangular",
                       "Circular"    = "circular"),
          selected = "rectangular"
        ),
        checkboxInput("tree_show_labels",  "Show tip labels",   value = TRUE),
        checkboxInput("tree_show_support", "Show node support", value = FALSE),
        # Tip color mode — replaces single checkbox
        radioButtons(
          "tree_color_by",
          "Color tips by:",
          choices  = c("None"      = "none",
                       "Country"   = "country",
                       "Haplotype" = "haplotype"),
          selected = "country"
        ),
        hr(class = "sidebar-divider"),
        tags$p(class = "help-text",
               icon("circle-info"),
               " Tip colors correspond to selected grouping.")
      ),
      # Right: tree plot
      card(
        full_screen = TRUE,
        card_header("RAxML Maximum Likelihood Tree — COI"),
        withSpinner(
          plotlyOutput("phylo_tree", height = "700px"),
          color = "#1a4f72"
        )
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
          tags$h4("About the AnStep Atlas"),
          tags$p(
            "This dashboard visualizes COI haplotype diversity and geographic distribution of ",
            tags$em("Anopheles stephensi"), " and related vectors across native and invasive ranges."
          ),
          tags$h5("Data"),
          tags$p(
            "Sequences were aligned with MAFFT and haplotypes were inferred from previous findings. ",
            "Phylogenetic inference was performed with RAxML. ", "Scripts available here: ", tags$a(href= "https://github.com/sayalg/AnStep_Atlas", "https://github.com/sayalg/AnStep_Atlas")
          ),
          tags$h5("Contact"),
          # TODO: update with your actual contact info
          tags$p("For questions or data requests, contact: ",
                 tags$a(href = "mailto:Tamar_Carter@baylor.edu", "Tamar_Carter@baylor.edu")),
          tags$h5("Citation"),
          tags$p(tags$em("TODO: Add citation once published.")),
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

  # ---- Palette (computed once) ----
  haplotype_palette <- create_haplotype_palette(df_haplo, category_col = "Haplotype")

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

  # Base map with static haplotype legend — renders once
  output$map <- renderLeaflet({
    create_base_map() %>%
      add_country_status_layer(
        data        = invasive_status,
        country_col = "COUNTRY_NAME",
        status_col  = "INVASIVE_STATUS"
      ) %>%
      add_railway_layer() %>%
      addLegend(
        position = "bottomleft",
        colors   = haplotype_palette,
        labels   = names(haplotype_palette),
        title    = "Haplotypes",
        opacity  = 1,
        layerId  = "haplotype_legend"
      )
  })

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

  # Update shipping routes layer when year filter changes
  filtered_shipping_sf <- reactive({
    if (is.null(shipping_routes_by_year) || length(shipping_years) == 0) return(NULL)

    target_years <- switch(input$time_mode,
      "exact"      = shipping_years[shipping_years == input$year_filter],
      "cumulative" = shipping_years[shipping_years <= input$year_filter],
      "range"      = shipping_years[shipping_years >= input$year_range[1] &
                                    shipping_years <= input$year_range[2]],
      "all"        = shipping_years
    )

    keys <- intersect(as.character(target_years), names(shipping_routes_by_year))
    if (length(keys) == 0) return(NULL)

    do.call(rbind, shipping_routes_by_year[keys])
  })

  observe({
    add_shipping_layer(leafletProxy("map"), filtered_shipping_sf())
  })

  # ==============================================================================
  # TAB 2 — PHYLOGENETIC TREE
  # ==============================================================================

  output$phylo_tree <- renderPlotly({
    tree_path <- "../RAXML/COI_subsampled20_per_country_year.aligned.raxml.bestTree"
    
    tree <- ape::read.tree(tree_path)
    
    layout_name <- switch(input$tree_layout,
                          "rectangular" = "rectangular",
                          "circular"    = "circular"
    )
    
    # Get tree coordinates directly from ggtree
    tree_df <- ggtree(tree, layout = layout_name)$data
    
    # Join metadata
    meta <- df_haplo %>% select(FASTA.ID, Country, Haplotype) %>% distinct()
    tree_df <- tree_df %>% left_join(meta, by = c("label" = "FASTA.ID"))
    
    tips  <- tree_df %>% filter(isTip == TRUE)
    nodes <- tree_df %>% filter(isTip == FALSE)
    
    # Color column
    color_by <- if (input$tree_color_by == "country")   tips$Country
    else if (input$tree_color_by == "haplotype") tips$Haplotype
    else NULL
    
    # Build edge traces (branches)
    edge_traces <- lapply(1:nrow(tree$edge), function(i) {
      parent <- tree_df[tree_df$node == tree$edge[i, 1], ]
      child  <- tree_df[tree_df$node == tree$edge[i, 2], ]
      list(
        x    = c(parent$x, child$x, NA),
        y    = c(parent$y, child$y, NA),
        type = "scatter", mode = "lines",
        line = list(color = "#2e86ab", width = 0.8),
        hoverinfo = "none", showlegend = FALSE
      )
    })
    
    pl <- plot_ly()
    
    # Add branches
    for (tr in edge_traces) {
      pl <- pl %>% add_trace(
        x = tr$x, y = tr$y,
        type = "scatter", mode = "lines",
        line = list(color = "#2e86ab", width = 0.8),
        hoverinfo = "none", showlegend = FALSE
      )
    }
    
    # Add tip points
    if (!is.null(color_by)) {
      tips$color_group <- color_by
      for (grp in unique(color_by)) {
        grp_tips <- tips %>% filter(color_group == grp)
        grp_hover <- paste0("ID: ", grp_tips$label, "<br>", input$tree_color_by, ": ", grp)
        pl <- pl %>% add_trace(
          x    = grp_tips$x, y = grp_tips$y,
          type = "scatter", mode = "markers",
          marker = list(size = 6),
          name   = grp,
          text   = grp_hover, hoverinfo = "text"
        )
      }
    } else {
      pl <- pl %>% add_trace(
        x = tips$x, y = tips$y,
        type = "scatter", mode = "markers",
        marker = list(size = 6, color = "#1a4f72"),
        text = paste0("ID: ", tips$label), hoverinfo = "text",
        showlegend = FALSE
      )
    }
    
    # Add tip labels
    if (input$tree_show_labels) {
      pl <- pl %>% add_annotations(
        x = tips$x, y = tips$y,
        text      = tips$label,
        showarrow = FALSE,
        xanchor   = "left",
        xshift    = 6,
        font      = list(size = 8, color = "#333")
      )
    }
    
    # Add node support labels
    if (input$tree_show_support) {
      pl <- pl %>% add_annotations(
        x = nodes$x, y = nodes$y,
        text      = nodes$label,
        showarrow = FALSE,
        xanchor   = "right",
        font      = list(size = 7, color = "#c0392b")
      )
    }
    
    pl %>% layout(
      xaxis      = list(title = "", showgrid = FALSE, zeroline = FALSE),
      yaxis      = list(title = "", showgrid = FALSE, zeroline = FALSE),
      dragmode   = "zoom",
      hoverlabel = list(bgcolor = "white"),
      plot_bgcolor  = "white",
      paper_bgcolor = "white"
    ) %>%
      config(
        scrollZoom             = TRUE,
        displaylogo            = FALSE,
        modeBarButtonsToRemove = c("lasso2d", "select2d")
      )
  })

  # ==============================================================================
  # TAB 3, PANEL 1 — STRAINHUB NETWORK
  # ==============================================================================

  output$haplo_network <- renderVisNetwork({

    # --- Minimal placeholder network ---
    treedata <- ape::read.tree("../RAXML/COI_subsampled20_per_country_year.aligned.raxml.bestTree")
    metadata <- readr::read_csv("../Input_Files/FinalCOI_metadata_w_coor.csv", col_names = TRUE)
    
    ## Adjust column names
    names(metadata)[names(metadata) == "FASTA.ID"] <- "Accession"
    names(metadata)[names(metadata) == "Accession.ID"] <- "Unique"
    
    ## Prune Tree
    taxa_without_data <- setdiff(treedata$tip.label, metadata$Accession)
    pruned_tree <- tidytree::drop.tip(treedata, taxa_without_data)
    
    ## Create Country with Haplotype column
    metadata$haplo_concat <- paste(metadata$Country, metadata$Haplotype, sep = "_")
    
    ## Make the Transmission Network
    graph <- make_transnet(pruned_tree,
                           metadata,
                           columnSelection = "Country",
                           # columnSelection = "haplo_concat",
                           centralityMetric = 3, #Betweenness Centrality
                           treeType = "parsimonious")
    
    ###############################################################################
    # Customize Strainhub Network Graph
    ###############################################################################
    
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
             color = "black",
             width = case_when(value == 1 ~ 1, value == 2 ~ 5, value == 3 ~ 20),
             value = NULL,
             length = 500)
    
    visNetwork(nodes, edges,
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

    metadata_cn <- readr::read_csv("../Input_Files/FinalCOI_metadata_w_coor.csv",
                                   col_names = TRUE, show_col_types = FALSE)

    # Clean metadata
    metadata_clean_cn <- metadata_cn %>%
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

    # Edges
    edges_idx    <- which(country_matrix > 0, arr.ind = TRUE)
    edges_idx    <- edges_idx[edges_idx[, 1] < edges_idx[, 2], ]
    links_country <- data.frame(
      from  = edges_idx[, 1],
      to    = edges_idx[, 2],
      value = mapply(function(f, t) country_matrix[f, t],
                     edges_idx[, 1], edges_idx[, 2])
    )

    visNetwork(nodes_country, links_country,
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
