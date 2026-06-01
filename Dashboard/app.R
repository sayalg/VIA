# ---- Load Required Packages ----
library(shiny)
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

# ---- Configuration ----
MAPBOX_TOKEN      <- "pk.eyJ1Ijoic2F5YWxnIiwiYSI6ImNtaWdha3Y4eDA1YmczZXEybjZvZjE0YTQifQ.T_HqxBxbmcdViI2L0LIzgQ"
MAPBOX_STYLE      <- "https://api.mapbox.com/styles/v1/sayalg/cmkr2q0c3000q01s87umrehau/tiles/256/{z}/{x}/{y}@2x?access_token={accessToken}"
CARTO_TILES       <- "https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png"
INVALID_LOCATIONS <- c("Laboratory")

# ================================================================================
# DATA LOADING AND PREPROCESSING
# ================================================================================

# ---- Load Data ----
samples_df <- read.csv("../Input_Files/FinalCOI_metadata.csv", sep = ",", header = TRUE)
samples_df[samples_df == ""] <- NA

invasive_status <- read_excel("../Input_Files/MTM_INVASIVE_VECTOR_SPECIES_20251205.xlsx", sheet = "Data")
invasive_status[invasive_status == ""] <- NA

# ---- Clean Data ----
samples_df <- samples_df %>%
  filter(
    !str_detect(City.Town, paste(INVALID_LOCATIONS, collapse = "|")) | is.na(City.Town),
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
# write.csv(df_geo, "../Input_Files/FinalCOI_metadata_w_coor.csv", row.names = FALSE)
df_geo <- read.csv("../Input_Files/FinalCOI_metadata_w_coor.csv", sep = ",", header = TRUE)

# Filter to valid haplotypes
df_haplo <- df_geo %>% filter(!is.na(Haplotype))

# ---- Cache World Polygons ----
# Loaded once at startup and passed into add_country_status_layer to avoid
# redundant downloads/reads on every map render or filter change.
world_polygons <- ne_countries(returnclass = "sf") %>%
  mutate(name_clean = tolower(admin))

# ================================================================================
# COLOR UTILITIES
# ================================================================================

# ---- High-Contrast Palette Generator ----
# Uses a hand-curated set of perceptually distinct colors so that adjacent
# haplotypes (e.g. Hap1/Hap2) are never visually similar. Falls back to
# interpolation for datasets with more than 24 haplotypes.
create_high_contrast_palette <- function(n) {
  base_colors <- c(
    "#E41A1C", "#4DAF4A", "#377EB8", "#FF7F00", "#984EA3",
    "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
    "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
    "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666"
  )
  if (n <= length(base_colors)) {
    return(base_colors[seq_len(n)])
  }
  colorRampPalette(base_colors)(n)
}

# ---- Build Named Haplotype Palette ----
create_haplotype_palette <- function(data, category_col = "Haplotype") {
  all_haplotypes <- data %>%
    filter(!is.na(.data[[category_col]])) %>%
    pull(.data[[category_col]]) %>%
    unique()

  # Sort numerically so Hap1, Hap2, ..., Hap10 not Hap1, Hap10, Hap2
  numeric_parts <- as.numeric(str_extract(all_haplotypes, "\\d+"))
  if (!all(is.na(numeric_parts))) {
    all_haplotypes <- all_haplotypes[order(numeric_parts, na.last = TRUE)]
  }

  colors <- create_high_contrast_palette(length(all_haplotypes))
  setNames(as.character(colors), all_haplotypes)
}

# ================================================================================
# MAP LAYER FUNCTIONS
# ================================================================================

# ---- Create Base Map ----
create_base_map <- function(tile_size = 256) {
  leaflet() %>%
    addMapPane("countryStatus", zIndex = 200) %>%
    addMapPane("charts",        zIndex = 400) %>%
    addProviderTiles(providers$CartoDB.Positron,   group = "Base") %>%
    addProviderTiles(providers$Esri.WorldImagery,  group = "Satellite") %>%
    addTiles(
      urlTemplate = MAPBOX_STYLE,
      options     = tileOptions(accessToken = MAPBOX_TOKEN, tileSize = tile_size),
      group       = "Land Use"
    ) %>%
    addLayersControl(
      baseGroups    = c("Base", "Satellite"),
      overlayGroups = c("Land Use", "Country Status"),
      options       = layersControlOptions(collapsed = FALSE)
    )
}

# ---- Add Pie Chart Markers ----
add_pie_markers <- function(map, data,
                            lat_col      = "Latitude",
                            lon_col      = "Longitude",
                            category_col = "Haplotype",
                            city_col     = "City.Town",
                            state_col    = "County.State",
                            country_col  = "Country",
                            radius       = 20,
                            add_legend   = TRUE,
                            color_palette = NULL) {

  # --- Location label: City → State → Country fallback ---
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

  # --- Aggregate haplotype counts by location ---
  df_pie <- data %>%
    filter(!is.na(.data[[lat_col]]), !is.na(.data[[lon_col]])) %>%
    count(.data[[lat_col]], .data[[lon_col]], .data[[category_col]]) %>%
    pivot_wider(
      names_from  = all_of(category_col),
      values_from = n,
      values_fill = 0
    ) %>%
    left_join(city_lookup, by = c(lat_col, lon_col))

  coords <- df_pie %>% select(all_of(c(lat_col, lon_col)))
  pies   <- df_pie %>% select(-all_of(c(lat_col, lon_col, "city")))

  # Sort haplotype columns numerically
  numeric_parts <- as.numeric(str_extract(colnames(pies), "\\d+"))
  if (!all(is.na(numeric_parts))) {
    pies <- pies[, order(numeric_parts, na.last = TRUE)]
  }

  # --- Resolve colors ---
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

  # --- Build per-location popup HTML ---
  layer_ids  <- make.unique(df_pie$city)
  color_map  <- setNames(colors, colnames(pies))
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
      "<span style='display:inline-block;width:10px;height:10px;",
      "background-color:", color_map[names(present)], ";",
      "border-radius:50%;margin-right:5px'></span>",
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
      rows_html,
      "</table></div>"
    )
  }

  # Convert named list → character vector in layer_id order for addMinicharts
  popup_html <- unlist(popup_list[layer_ids])

  # --- Add minichart pies to map ---
  # Passing `popup` directly to addMinicharts overrides its built-in popup
  # (which would show all haplotypes including zeros) with our pre-filtered HTML.
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

  # --- Optional static legend ---
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

# ---- Add Country Status Layer ----
# Accepts a pre-loaded `world` sf object to avoid redundant ne_countries() calls.
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

# ---- Fix Chart Pane Z-Index ----
# Moves minichart containers into the dedicated "charts" pane so they always
# render above country-status polygons. Click/popup logic lives in
# add_pie_markers via its own onRender call.
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

# ================================================================================
# SHINY APP
# ================================================================================

run_map_with_year_filter <- function(data          = df_haplo,
                                     invasive_data = invasive_status,
                                     year_col      = "Year") {

  if (!year_col %in% colnames(data)) {
    stop(paste("Column", year_col, "not found in data"))
  }

  available_years <- data %>%
    filter(!is.na(.data[[year_col]])) %>%
    pull(.data[[year_col]]) %>%
    unique() %>%
    sort()

  if (length(available_years) == 0) {
    stop("No valid years found in the data")
  }

  # --- UI ---
  ui <- fluidPage(
    titlePanel("AnStep Atlas"),

    sidebarLayout(
      sidebarPanel(
        width = 3,

        radioButtons(
          "time_mode", "Filter Mode:",
          choices  = c("Exact Year" = "exact",
                       "Up to Year" = "cumulative",
                       "All Years"  = "all"),
          selected = "all"
        ),

        # Slider is hidden when "All Years" is selected
        conditionalPanel(
          condition = "input.time_mode !== 'all'",
          sliderInput(
            "year_filter", "Select Year:",
            min     = min(available_years),
            max     = max(available_years),
            value   = max(available_years),
            step    = 1,
            sep     = "",
            animate = animationOptions(interval = 1000, loop = FALSE)
          )
        ),

        hr(),
        textOutput("sample_count"),
        uiOutput("no_data_warning")
      ),

      mainPanel(
        width = 9,
        withSpinner(leafletOutput("map", height = "800px"))
      )
    )
  )

  # --- Server ---
  server <- function(input, output, session) {

    # Computed once at startup — the palette covers all haplotypes in the full
    # dataset so colors remain consistent regardless of which years are shown.
    haplotype_palette <- create_haplotype_palette(data, category_col = "Haplotype")

    # Filter data according to the selected time mode
    filtered_data <- reactive({
      switch(input$time_mode,
        "exact"      = data %>% filter(.data[[year_col]] == input$year_filter),
        "cumulative" = data %>% filter(.data[[year_col]] <= input$year_filter),
        "all"        = data
      )
    })

    # Sample count summary line
    output$sample_count <- renderText({
      n <- nrow(filtered_data())
      s <- if (n != 1) "s" else ""
      switch(input$time_mode,
        "exact"      = paste0(n, " sample", s, " from ", input$year_filter),
        "cumulative" = paste0(n, " sample", s, " up to ", input$year_filter),
        "all"        = paste0(n, " sample", s, " across all years")
      )
    })

    # Warning banner shown when a filter produces no results
    output$no_data_warning <- renderUI({
      if (nrow(filtered_data()) == 0 && input$time_mode != "all") {
        tags$div(
          style = paste0(
            "margin-top:10px;padding:8px 12px;",
            "background:#fff3cd;border:1px solid #ffc107;",
            "border-radius:4px;color:#856404;font-size:13px"
          ),
          icon("exclamation-triangle"),
          " No samples found for the selected year(s)."
        )
      }
    })

    # Base map with static haplotype legend (rendered once; only pie charts
    # are updated as the filter changes)
    output$map <- renderLeaflet({
      create_base_map() %>%
        add_country_status_layer(
          data        = invasive_data,
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

    # Update pie charts whenever the filter changes
    observe({
      req(filtered_data())
      fd    <- filtered_data()
      proxy <- leafletProxy("map") %>% clearMinicharts()

      if (nrow(fd) > 0) {
        proxy <- add_pie_markers(
          map          = proxy,
          data         = fd,
          lat_col      = "Latitude",
          lon_col      = "Longitude",
          category_col = "Haplotype",
          city_col     = "City.Town",
          state_col    = "County.State",
          country_col  = "Country",
          radius       = 30,
          add_legend   = FALSE,    # Legend is static on the base map
          color_palette = haplotype_palette
        )
      }

      add_chart_interactions(proxy)
    })
  }

  shinyApp(ui = ui, server = server)
}

run_map_with_year_filter()
