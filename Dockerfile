# VIA: Vector Invasion Atlas, as a container.
#
#   docker build -t via .
#   docker run --rm -p 3838:3838 via          # then open http://localhost:3838
#
# Optional CARTO basemap (the app uses a keyless Esri basemap without it):
#   docker run --rm -p 3838:3838 -e CARTO_API_KEY=your_key via
#
# The image holds only the app: Dashboard/ with the data the build already
# produced (Dashboard/data, Dashboard/www). To change the data, rerun
# Tools/build_app_cache.R on the host and rebuild the image.

# R 4.5.2 matches Dashboard/manifest.json. rocker/r-ver installs packages from a
# dated Posit Package Manager snapshot, so rebuilding later gives the same
# package versions, as Linux binaries (no compiling).
FROM rocker/r-ver:4.5.2

# System libraries: sf needs GDAL/GEOS/PROJ/udunits; httr, xml2 and the web
# stack need curl, openssl and libxml2.
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl \
        libabsl-dev \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libfribidi-dev \
        libgdal-dev \
        libgeos-dev \
        libharfbuzz-dev \
        libjpeg-dev \
        libpng-dev \
        libproj-dev \
        libsqlite3-dev \
        libssl-dev \
        libtiff-dev \
        libudunits2-dev \
        libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Only the packages the running app loads (the data build's extras are not
# needed here). Dependencies come along.
RUN install2.r --error --skipinstalled \
        shiny bslib shinycssloaders leaflet leaflet.minicharts \
        dplyr tidyr sf visNetwork stringr httr xml2 \
    && rm -rf /tmp/downloaded_packages

# The app reads "data/" and "www/" relative to its own folder.
WORKDIR /app
COPY Dashboard/ /app/

# Don't run as root.
RUN useradd --system --create-home shiny
USER shiny

EXPOSE 3838
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -fs http://localhost:3838/ >/dev/null || exit 1

CMD ["R", "-q", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838)"]
