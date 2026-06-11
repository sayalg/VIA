# ================================================================================
# prepare_shipping_routes.R
# Run ONCE offline.
#
# Filters routes to those originating from the Horn of Africa, Middle East,
# or India. Routes are drawn as multi-segment LINESTRING geometries that pass
# through maritime waypoints, keeping lines over water.
#
# Input:  ../Input_Files/shipping_routes.csv  (tab-separated)
# Output: ../Input_Files/shipping_routes.rds
#         Named list: list("2006" = sf, "2007" = sf, ...)
# ================================================================================

library(sf)
library(tidyverse)
library(rnaturalearth)
library(reticulate)

OUT_RDS <- "../Input_Files/shipping_routes.rds"

# ================================================================================
# REGION FILTER
# Countries of origin to keep (Horn of Africa + Middle East + India)
# ================================================================================
origin_countries <- c(
  # Horn of Africa
  "Somalia", "Ethiopia", "Eritrea", "Djibouti", "Kenya",
  # Middle East
  "Yemen", "Oman", "United Arab Emirates", "Saudi Arabia", "Qatar",
  "Bahrain", "Kuwait", "Iraq", "Iran", "Jordan", "Israel", "Lebanon",
  "Syria", "Egypt",
  # India
  "India"
)

# ================================================================================
# MARITIME WAYPOINTS
# Named lon/lat points at key straits and open-ocean nodes.
# Routes are stitched together through relevant waypoints so they stay over water.
# ================================================================================
wp <- list(
  bab_el_mandeb  = c(43.4,  12.6),   # Strait between Yemen and Djibouti
  gulf_of_aden   = c(48.0,  12.0),   # Mid Gulf of Aden
  suez_s         = c(32.6,  29.9),   # Southern Suez Canal entrance
  suez_n         = c(32.3,  31.2),   # Northern Suez Canal entrance
  red_sea_mid    = c(38.0,  20.0),   # Mid Red Sea
  hormuz         = c(56.5,  26.2),   # Strait of Hormuz
  gulf_oman      = c(59.0,  23.5),   # Gulf of Oman open water
  arabian_sea_w  = c(58.0,  15.0),   # Western Arabian Sea
  arabian_sea_e  = c(66.0,  15.0),   # Eastern Arabian Sea
  india_sw       = c(74.0,   9.0),   # SW India coast
  india_se       = c(80.5,   8.5),   # SE India / Palk Strait area
  india_nw       = c(70.0,  20.0),   # NW India (Gujarat coast)
  malacca        = c(103.8,  1.3),   # Strait of Malacca
  med_e          = c(28.0,  34.5),   # Eastern Mediterranean
  med_w          = c(5.0,   37.0),   # Western Mediterranean
  gibraltar      = c(-5.4,  35.9),   # Strait of Gibraltar
  cape_good_hope = c(18.4, -34.4),   # Cape of Good Hope
  e_africa_coast = c(45.0,  -5.0),   # East African coast
  indian_ocean_s = c(70.0, -20.0),   # Southern Indian Ocean
  se_asia        = c(108.0,  5.0),   # South China Sea entry
  pacific_w      = c(130.0, 20.0)    # Western Pacific
)

# Helper: return lon/lat matrix for a named waypoint
W <- function(...) {
  pts <- list(...)
  do.call(rbind, lapply(pts, function(p) wp[[p]]))
}

# ================================================================================
# ROUTING LOGIC
# Assign waypoints based on origin and destination broad region.
# Waypoints are inserted between origin and destination centroid.
# ================================================================================

broad_region <- function(lon, lat) {
  # Classify a lon/lat into a broad ocean-routing region
  if (lat > 25  & lon > 25  & lon < 60)  return("gulf")          # Persian Gulf / Red Sea
  if (lat > 5   & lon > 60  & lon < 85)  return("india_w")       # W India / Arabian Sea
  if (lat > 5   & lon > 85  & lon < 100) return("india_e")       # E India / Bay of Bengal
  if (lon > 100 & lat < 25)              return("se_asia")        # SE Asia
  if (lon > 100 & lat > 25)              return("east_asia")      # E Asia
  if (lat > 30  & lon > 10  & lon < 40)  return("med_e")         # Eastern Med
  if (lat > 30  & lon < 10)             return("atlantic")        # Atlantic / N Europe
  if (lat < 0   & lon > 10  & lon < 55)  return("e_africa_s")    # S East Africa
  if (lat < -20 & lon > 55)             return("indian_ocean_s")  # S Indian Ocean
  return("other")
}

# Returns a list of intermediate waypoint matrices for a given origin → dest region pair
get_waypoints <- function(o_region, d_region, o_lon, o_lat, d_lon, d_lat) {

  # Routes from Gulf / Red Sea origins
  if (o_region == "gulf") {
    through_hormuz <- o_lon > 48   # Gulf countries go through Hormuz; Red Sea stays in Red Sea

    if (d_region == "india_w")     return(if(through_hormuz) W("hormuz","gulf_oman","arabian_sea_w") else W("bab_el_mandeb","arabian_sea_w"))
    if (d_region == "india_e")     return(if(through_hormuz) W("hormuz","gulf_oman","arabian_sea_e","india_sw") else W("bab_el_mandeb","arabian_sea_e","india_sw"))
    if (d_region == "se_asia")     return(if(through_hormuz) W("hormuz","gulf_oman","arabian_sea_e","india_sw","malacca") else W("bab_el_mandeb","arabian_sea_e","india_sw","malacca"))
    if (d_region == "east_asia")   return(if(through_hormuz) W("hormuz","gulf_oman","arabian_sea_e","india_sw","malacca","se_asia","pacific_w") else W("bab_el_mandeb","arabian_sea_e","india_sw","malacca","se_asia"))
    if (d_region == "med_e")       return(W("red_sea_mid","suez_s","suez_n","med_e"))
    if (d_region == "atlantic")    return(W("red_sea_mid","suez_s","suez_n","med_e","med_w","gibraltar"))
    if (d_region == "e_africa_s")  return(W("bab_el_mandeb","e_africa_coast"))
    if (d_region == "indian_ocean_s") return(W("bab_el_mandeb","arabian_sea_w","indian_ocean_s"))
  }

  # Routes from India (west coast)
  if (o_region == "india_w") {
    if (d_region == "gulf")        return(W("arabian_sea_w","hormuz"))
    if (d_region == "se_asia")     return(W("india_sw","malacca"))
    if (d_region == "east_asia")   return(W("india_sw","malacca","se_asia"))
    if (d_region == "med_e")       return(W("arabian_sea_w","bab_el_mandeb","red_sea_mid","suez_s","suez_n","med_e"))
    if (d_region == "atlantic")    return(W("arabian_sea_w","bab_el_mandeb","red_sea_mid","suez_s","suez_n","med_w","gibraltar"))
    if (d_region == "e_africa_s")  return(W("arabian_sea_w","e_africa_coast"))
    if (d_region == "indian_ocean_s") return(W("arabian_sea_w","indian_ocean_s"))
  }

  # Routes from India (east coast)
  if (o_region == "india_e") {
    if (d_region == "gulf")        return(W("india_sw","arabian_sea_w","hormuz"))
    if (d_region == "se_asia")     return(W("malacca"))
    if (d_region == "east_asia")   return(W("malacca","se_asia"))
    if (d_region == "med_e")       return(W("india_sw","arabian_sea_w","bab_el_mandeb","red_sea_mid","suez_s","suez_n","med_e"))
    if (d_region == "atlantic")    return(W("india_sw","arabian_sea_w","bab_el_mandeb","red_sea_mid","suez_s","suez_n","med_w","gibraltar"))
    if (d_region == "e_africa_s")  return(W("india_sw","e_africa_coast"))
  }

  return(NULL)   # fallback: straight line (same region or unhandled)
}

# ================================================================================
# BUILD ROUTED LINESTRING via searoute-py
# ================================================================================
sr <- import("searoute")

make_routed_line <- function(from_lon, from_lat, to_lon, to_lat) {
  route <- tryCatch(
    sr$searoute(list(from_lon, from_lat), list(to_lon, to_lat)),
    error = function(e) NULL
  )
  
  # Fallback to straight line if searoute fails
  if (is.null(route)) {
    coords <- rbind(c(from_lon, from_lat), c(to_lon, to_lat))
    return(st_linestring(coords))
  }
  
  # Extract [lon, lat] coordinate pairs from GeoJSON output
  coords <- do.call(rbind, lapply(route$geometry$coordinates, function(p) c(p[[1]], p[[2]])))
  st_linestring(coords)
}

# ================================================================================
# 1. Country centroids
# ================================================================================
message("Building country centroids...")
sf::sf_use_s2(FALSE)
centroids <- ne_countries(returnclass = "sf") %>%
  st_centroid() %>%
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2]) %>%
  st_drop_geometry() %>%
  select(name = admin, lon, lat)
sf::sf_use_s2(TRUE)
message(sprintf("  %d centroids loaded.", nrow(centroids)))

# ================================================================================
# 2. Load CSV, filter to origin region, join coordinates
# ================================================================================
message("Reading and filtering CSV...")
raw <- read.csv("US_LSBCI.csv",
                sep = ",", header = TRUE, check.names = FALSE)

routes_df <- raw %>%
  rename(from = `Economy Label`, to = `Partner Label`, quarter = `Quarter Label`) %>%
  mutate(year = as.integer(str_extract(quarter, "\\d{4}"))) %>%
  filter(!is.na(year), from != to) %>%
  filter(from %in% origin_countries) %>%        # ← region filter
  distinct(from, to, year) %>%
  inner_join(centroids, by = c("from" = "name")) %>%
  rename(from_lon = lon, from_lat = lat) %>%
  inner_join(centroids, by = c("to" = "name")) %>%
  rename(to_lon = lon, to_lat = lat)

message(sprintf("  %d unique from-to-year combinations after filtering.", nrow(routes_df)))

# ================================================================================
# 3. Build routed sf LINESTRINGs per year
# ================================================================================
message("Building routed geometries by year...")

build_sf_for_year <- function(df) {
  lines <- lapply(seq_len(nrow(df)), function(i) {
    make_routed_line(df$from_lon[i], df$from_lat[i],
                     df$to_lon[i],   df$to_lat[i])
  })
  st_sf(from     = df$from,
        to       = df$to,
        geometry = st_sfc(lines, crs = 4326))
}

years <- sort(unique(routes_df$year))
shipping_by_year <- setNames(
  lapply(years, function(y) {
    sf_y <- build_sf_for_year(filter(routes_df, year == y))
    message(sprintf("  %d: %d routes", y, nrow(sf_y)))
    sf_y
  }),
  as.character(years)
)

# ================================================================================
# 4. Save
# ================================================================================
message(sprintf("Saving to %s ...", OUT_RDS))
dir.create(dirname(OUT_RDS), showWarnings = FALSE, recursive = TRUE)
saveRDS(shipping_by_year, OUT_RDS)
message("Done.")
