# ================================================================================
# prepare_africa_railways.R
# Run this ONCE offline to download all African country railway shapefiles from
# Geofabrik and merge them into a single cached RDS file.
#
# Output: ../Dashboard/data/africa_railways.rds  (sf object, CRS = WGS84)
#         This is read directly by Dashboard/app.R at startup.
#
# After running, do NOT re-run unless you want to refresh the OSM data.
# ================================================================================

library(httr)
library(sf)
library(tidyverse)

# ---- Config ----
RAILWAY_SHP <- "gis_osm_railways_free_1.shp"
OUT_RDS     <- "../Dashboard/data/africa_railways.rds"
TMP_DIR     <- tempdir()

# ---- 1. Build URLs from hardcoded country list ----
# Geofabrik uses lowercase hyphenated slugs
african_countries <- c(
  "algeria", "angola", "benin", "botswana", "burkina-faso", "burundi",
  "cameroon", "canary-islands", "cape-verde", "central-african-republic",
  "chad", "comores", "congo-brazzaville", "congo-democratic-republic",
  "djibouti", "egypt", "equatorial-guinea", "eritrea", "ethiopia",
  "gabon", "ghana", "guinea", "guinea-bissau", "ivory-coast",
  "kenya", "lesotho", "liberia", "libya", "madagascar", "malawi",
  "mali", "mauritania", "mauritius", "morocco", "mozambique", "namibia",
  "niger", "nigeria", "rwanda", "saint-helena-ascension-and-tristan-da-cunha",
  "sao-tome-and-principe", "senegal-and-gambia", "seychelles",
  "sierra-leone", "somalia", "south-africa", "south-sudan", "sudan",
  "swaziland", "tanzania", "togo", "tunisia", "uganda",
  "western-sahara", "zambia", "zimbabwe"
)

country_zips <- paste0(
  "https://download.geofabrik.de/africa/",
  african_countries,
  "-latest-free.shp.zip"
)

# ---- 2. Validate URLs with HEAD requests ----
message(sprintf("Validating %d URLs...", length(country_zips)))
valid <- vapply(country_zips, function(url) {
  tryCatch({
    code <- httr::status_code(httr::HEAD(url))
    if (code != 200) message(sprintf("  SKIP (HTTP %d): %s", code, basename(url)))
    code == 200
  }, error = function(e) {
    message(sprintf("  SKIP (error): %s", basename(url)))
    FALSE
  })
}, logical(1))

country_zips <- country_zips[valid]
message(sprintf("%d / %d URLs valid.", sum(valid), length(valid)))
if (length(country_zips) == 0) stop("No valid URLs found.")

# ---- 3. Download, extract, read railway layer ----
railway_list <- list()

for (url in country_zips) {
  country_name <- sub(".*africa/(.+)-latest-free\\.shp\\.zip", "\\1", url)
  zip_path     <- file.path(TMP_DIR, basename(url))
  extract_dir  <- file.path(TMP_DIR, country_name)

  message(sprintf("Downloading %s ...", country_name))
  tryCatch({
    download.file(url, destfile = zip_path, mode = "wb", quiet = TRUE)
    dir.create(extract_dir, showWarnings = FALSE)

    shp_files <- paste0(sub("\\.shp$", "", RAILWAY_SHP), c(".shp", ".dbf", ".shx", ".prj"))
    unzip(zip_path, files = shp_files, exdir = extract_dir, overwrite = TRUE)

    shp_path <- file.path(extract_dir, RAILWAY_SHP)
    if (file.exists(shp_path)) {
      sf_obj <- st_read(shp_path, quiet = TRUE)
      sf_obj$country <- country_name
      railway_list[[country_name]] <- sf_obj
      message(sprintf("  -> %d features", nrow(sf_obj)))
    } else {
      message("  -> railway shapefile not found in zip, skipping.")
    }
  }, error = function(e) {
    message(sprintf("  ERROR: %s", conditionMessage(e)))
  })
}

# ---- 4. Merge, filter, and save ----
if (length(railway_list) == 0) stop("No railway data downloaded.")

message("Merging all countries...")
africa_railways <- bind_rows(railway_list) |>
  st_as_sf() |>
  st_transform(4326) |>
  filter(fclass %in% c("rail", "light_rail", "narrow_gauge", "preserved"))

message(sprintf("Total railway features after filtering: %d", nrow(africa_railways)))

dir.create(dirname(OUT_RDS), showWarnings = FALSE, recursive = TRUE)
saveRDS(africa_railways, OUT_RDS)
message(sprintf("Done. Saved to: %s", OUT_RDS))
