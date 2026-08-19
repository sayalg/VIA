# =============================================================================
# ANOPHELES STEPHENSI PUBMED PUBLICATION BROWSER — SHINY MODULE (v3)
# =============================================================================
# HOW TO INTEGRATE INTO YOUR EXISTING SHINY APP:
#
#   STEP 1 — Place this file (pubmed_module.R) in your app's folder
#   STEP 2 — At the top of your app.R:  source("pubmed_module.R")
#   STEP 3 — In your UI add a tab:      pubBrowserUI("pubs")
#   STEP 4 — In your server add a line: pubBrowserServer("pubs")
#   STEP 5 — Packages: install.packages(c("shiny","dplyr","stringr","httr","xml2"))
#
# WHAT CHANGED vs v2 (behaviour is identical, only efficiency improved):
#   * PubMed results are cached at the PROCESS level with a time-to-live, so the
#     20–40s network fetch happens once and is then shared by every user session
#     instead of firing on every session start. The ↻ Refresh button forces a
#     fresh fetch. This is the single biggest responsiveness win for the module.
#   * Keyword parsing in the card renderer is done once per row instead of
#     splitting the same string three times.
# =============================================================================

library(shiny)
library(dplyr)
library(stringr)
library(httr)
library(xml2)

# =============================================================================
# PROCESS-LEVEL CACHE  (shared across all sessions of this R process)
# =============================================================================

.pubmed_cache <- new.env(parent = emptyenv())
.pubmed_cache$data <- NULL          # data.frame of fetched papers
.pubmed_cache$time <- NULL          # POSIXct of last successful fetch

# Re-fetch from PubMed at most once per this window (seconds).
PUBMED_CACHE_TTL <- 6 * 60 * 60     # 6 hours

# =============================================================================
# PUBMED FETCH HELPERS
# =============================================================================

.fetch_pubmed_ids <- function(search_term, max_results = 500) {
  resp <- GET(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
    query = list(
      db      = "pubmed",
      term    = search_term,
      retmax  = max_results,
      retmode = "xml",
      sort    = "pub+date"
    )
  )
  if (http_error(resp)) return(character(0))
  doc <- read_xml(content(resp, as = "text", encoding = "UTF-8"))
  xml_text(xml_find_all(doc, "//Id"))
}

.fetch_pubmed_details <- function(ids) {
  if (length(ids) == 0) return(data.frame())
  batches  <- split(ids, ceiling(seq_along(ids) / 100))
  all_rows <- list()

  for (batch in batches) {
    resp <- GET(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
      query = list(
        db      = "pubmed",
        id      = paste(batch, collapse = ","),
        retmode = "xml"
      )
    )
    if (http_error(resp)) next
    doc      <- read_xml(content(resp, as = "text", encoding = "UTF-8"))
    articles <- xml_find_all(doc, "//PubmedArticle")

    # Build each batch's rows as a list, then bind once (avoids repeated
    # data.frame growth).
    batch_rows <- lapply(articles, function(art) {
      title <- xml_text(xml_find_first(art, ".//ArticleTitle"))
      title <- str_trim(str_remove_all(title, "\\[|\\]$"))

      author_nodes <- xml_find_all(art, ".//Author")
      authors <- vapply(author_nodes, function(a) {
        last <- xml_text(xml_find_first(a, "LastName"))
        init <- xml_text(xml_find_first(a, "Initials"))
        if (!is.na(last) && last != "") paste(last, init) else NA_character_
      }, character(1))
      authors <- paste(na.omit(authors), collapse = ", ")
      if (authors == "") authors <- "Unknown"

      year_node <- xml_find_first(art, ".//PubDate/Year")
      if (is.na(year_node)) year_node <- xml_find_first(art, ".//PubDate/MedlineDate")
      year <- str_extract(xml_text(year_node), "\\d{4}")
      if (is.na(year)) year <- "N/A"

      journal <- xml_text(xml_find_first(art, ".//Journal/Title"))
      if (is.na(journal) || journal == "")
        journal <- xml_text(xml_find_first(art, ".//ISOAbbreviation"))

      abstract <- paste(xml_text(xml_find_all(art, ".//AbstractText")), collapse = " ")
      abstract <- str_trim(abstract)
      if (abstract == "") abstract <- "No abstract available."

      doi_node <- xml_find_first(art, ".//ELocationID[@EIdType='doi']")
      doi      <- xml_text(doi_node)
      doi_url  <- if (!is.na(doi) && doi != "") {
        paste0("https://doi.org/", doi)
      } else {
        paste0("https://pubmed.ncbi.nlm.nih.gov/",
               xml_text(xml_find_first(art, ".//PMID")))
      }

      mesh_nodes <- xml_find_all(art, ".//MeshHeading/DescriptorName")
      kw_nodes   <- xml_find_all(art, ".//Keyword")
      keywords   <- unique(c(xml_text(mesh_nodes), xml_text(kw_nodes)))
      keywords   <- paste(keywords[seq_len(min(6, length(keywords)))], collapse = ", ")

      data.frame(
        title, authors, year, journal,
        doi = doi_url, abstract, keywords,
        stringsAsFactors = FALSE
      )
    })

    all_rows <- c(all_rows, batch_rows)
    Sys.sleep(0.34)   # respect NCBI ~3 req/s rate limit
  }

  if (length(all_rows) == 0) return(data.frame())
  bind_rows(all_rows)
}

# Fetch (or reuse cached) papers. Returns a data.frame, or NULL on failure.
.get_pubmed_papers <- function(force = FALSE) {
  fresh <- !is.null(.pubmed_cache$data) &&
    !is.null(.pubmed_cache$time) &&
    as.numeric(difftime(Sys.time(), .pubmed_cache$time, units = "secs")) < PUBMED_CACHE_TTL

  if (!force && fresh) return(.pubmed_cache$data)

  ids <- .fetch_pubmed_ids("Anopheles stephensi[Title]", 500)
  if (length(ids) == 0) return(NULL)
  df <- .fetch_pubmed_details(ids)
  if (nrow(df) == 0) return(NULL)

  .pubmed_cache$data <- df
  .pubmed_cache$time <- Sys.time()
  df
}

# =============================================================================
# MODULE CSS  (injected once when the UI loads)
# =============================================================================

.pubBrowserCSS <- function() {
  tags$style(HTML("
    /* Scoped under .pub-browser-wrap so styles don't bleed into your app */
    .pub-browser-wrap {
      background-color: #f4f6f9;
      min-height: 60vh;
      font-family: 'Source Sans 3', sans-serif;
      font-weight: 300;
    }

    /* ---- Header ---- */
    .pub-browser-wrap .pb-header {
      background: linear-gradient(135deg, #ffffff 0%, #eaf4f0 50%, #e8f0fa 100%);
      border-bottom: 1px solid #d8e4ed;
      padding: 36px 44px 28px;
      position: relative; overflow: hidden;
    }
    .pub-browser-wrap .pb-header::before {
      content: ''; position: absolute; top: -60px; right: -60px;
      width: 280px; height: 280px;
      background: radial-gradient(circle, rgba(46,138,94,0.10) 0%, transparent 70%);
      border-radius: 50%;
    }
    .pub-browser-wrap .pb-eyebrow {
      font-size: 11px; font-weight: 600; letter-spacing: 3px;
      text-transform: uppercase; color: #2e8a5e; margin-bottom: 8px;
    }
    .pub-browser-wrap .pb-title {
      font-size: 32px; font-weight: 700; color: #1a2535;
      line-height: 1.2; margin-bottom: 6px;
      font-family: Georgia, 'Times New Roman', serif;
    }
    .pub-browser-wrap .pb-title em { font-style: italic; color: #2e8a5e; }
    .pub-browser-wrap .pb-subtitle { font-size: 14px; color: #5a6f85; }

    /* ---- Controls ---- */
    .pub-browser-wrap .pb-controls {
      background: #ffffff;
      border-bottom: 1px solid #d8e4ed;
      padding: 12px 44px;
      display: flex; align-items: center; gap: 12px;
      box-shadow: 0 2px 10px rgba(0,0,0,0.06);
      flex-wrap: wrap;
    }
    .pub-browser-wrap .pb-search-icon { color: #2e8a5e; font-size: 17px; flex-shrink: 0; }
    .pub-browser-wrap .pb-controls .form-group { margin-bottom: 0 !important; flex: 1 1 220px; }

    .pub-browser-wrap input[id$='_search'] {
      width: 100% !important;
      background: #f4f6f9 !important; border: 1px solid #c8d8e8 !important;
      border-radius: 6px !important; color: #1a2535 !important;
      font-size: 14px !important; padding: 9px 14px !important;
      outline: none !important;
      transition: border-color 0.2s, box-shadow 0.2s !important;
    }
    .pub-browser-wrap input[id$='_search']:focus {
      border-color: #2e8a5e !important;
      box-shadow: 0 0 0 3px rgba(46,138,94,0.12) !important;
      background: #fff !important;
    }
    .pub-browser-wrap input[id$='_search']::placeholder { color: #9aaabb !important; }

    .pub-browser-wrap .pb-ctrl-group { display: flex; align-items: center; gap: 6px; flex-shrink: 0; }
    .pub-browser-wrap .pb-ctrl-label {
      font-size: 11px; color: #8a9daf; letter-spacing: 1px;
      text-transform: uppercase; white-space: nowrap;
    }
    .pub-browser-wrap .pb-ctrl-group .form-group { margin-bottom: 0 !important; }
    .pub-browser-wrap .pb-ctrl-group select,
    .pub-browser-wrap select[id$='_sort'] {
      background: #f4f6f9; border: 1px solid #c8d8e8; border-radius: 6px;
      color: #1a2535; font-size: 13px; padding: 7px 10px;
      cursor: pointer; outline: none;
    }
    .pub-browser-wrap .pb-ctrl-group select:focus,
    .pub-browser-wrap select[id$='_sort']:focus { border-color: #2e8a5e; }

    .pub-browser-wrap button[id$='_refresh'] {
      background: transparent; border: 1px solid #c8d8e8; border-radius: 6px;
      color: #2e8a5e; font-size: 11px; font-weight: 600; letter-spacing: 1px;
      text-transform: uppercase; padding: 7px 13px; cursor: pointer;
      transition: border-color 0.2s, background 0.2s; flex-shrink: 0;
      height: 34px;
    }
    .pub-browser-wrap button[id$='_refresh']:hover {
      background: rgba(46,138,94,0.08); border-color: #2e8a5e;
    }

    .pub-browser-wrap .pb-count {
      font-size: 12px; color: #8a9daf; white-space: nowrap; flex-shrink: 0;
    }
    .pub-browser-wrap .pb-count span { color: #2e8a5e; font-weight: 600; }

    /* ---- Loading / error ---- */
    .pub-browser-wrap .pb-loading,
    .pub-browser-wrap .pb-error {
      text-align: center; padding: 80px 20px;
    }
    .pub-browser-wrap .pb-spinner {
      width: 40px; height: 40px;
      border: 3px solid #d8e8f0; border-top-color: #2e8a5e;
      border-radius: 50%; animation: pb-spin 0.8s linear infinite;
      margin: 0 auto 18px;
    }
    @keyframes pb-spin { to { transform: rotate(360deg); } }
    .pub-browser-wrap .pb-loading p { font-size: 14px; color: #7a8f9f; }
    .pub-browser-wrap .pb-loading small { font-size: 12px; color: #9aaabb; display: block; margin-top: 5px; }
    .pub-browser-wrap .pb-error h3 { font-size: 18px; color: #c0392b; margin-bottom: 8px; }
    .pub-browser-wrap .pb-error p  { font-size: 13px; color: #8a5a55; }

    /* ---- Masonry grid ---- */
    .pub-browser-wrap .pb-grid {
      columns: 300px; column-gap: 16px; padding: 24px 40px 44px;
    }
    .pub-browser-wrap .pb-card-wrap { break-inside: avoid; margin-bottom: 16px; display: block; }

    /* ---- Card ---- */
    .pub-browser-wrap .pb-card {
      background: #ffffff; border: 1px solid #dce8f0; border-radius: 10px;
      padding: 20px; cursor: pointer;
      transition: transform 0.2s ease, border-color 0.2s ease, box-shadow 0.2s ease;
      display: flex; flex-direction: column; gap: 9px;
      position: relative; overflow: hidden; text-decoration: none;
    }
    .pub-browser-wrap .pb-card::before {
      content: ''; position: absolute; top: 0; left: 0;
      width: 3px; height: 100%;
      background: linear-gradient(180deg, #2e8a5e, #2e64a0);
      opacity: 0; transition: opacity 0.2s;
    }
    .pub-browser-wrap .pb-card:hover {
      transform: translateY(-3px); border-color: #a8c8e0;
      box-shadow: 0 8px 28px rgba(0,0,0,0.09);
    }
    .pub-browser-wrap .pb-card:hover::before { opacity: 1; }

    .pub-browser-wrap .pb-card-meta { display: flex; justify-content: space-between; align-items: center; }
    .pub-browser-wrap .pb-year {
      font-size: 10px; font-weight: 600; letter-spacing: 2px;
      color: #2e8a5e; text-transform: uppercase;
    }
    .pub-browser-wrap .pb-arrow { font-size: 13px; color: #c0ccd8; transition: color 0.2s; }
    .pub-browser-wrap .pb-card:hover .pb-arrow { color: #2e8a5e; }

    .pub-browser-wrap .pb-card-title {
      font-family: Georgia, 'Times New Roman', serif;
      font-size: 14px; font-weight: 600; color: #1a2535; line-height: 1.4;
    }
    .pub-browser-wrap .pb-authors { font-size: 12px; color: #6a7f95; font-style: italic; line-height: 1.4; }
    .pub-browser-wrap .pb-journal { font-size: 11px; font-weight: 600; color: #2e7a9a; letter-spacing: 0.4px; }
    .pub-browser-wrap .pb-abstract {
      font-size: 12px; color: #6a7f95; line-height: 1.5;
      display: -webkit-box; -webkit-line-clamp: 3;
      -webkit-box-orient: vertical; overflow: hidden;
    }
    .pub-browser-wrap .pb-tags { display: flex; flex-wrap: wrap; gap: 4px; margin-top: 2px; }
    .pub-browser-wrap .pb-tag {
      background: #eef4f8; border: 1px solid #c8dcea; border-radius: 20px;
      font-size: 9px; font-weight: 600; letter-spacing: 0.4px;
      color: #5a7a95; padding: 2px 7px; text-transform: uppercase;
    }

    /* ---- No results ---- */
    .pub-browser-wrap .pb-empty { text-align: center; padding: 60px 20px; }
    .pub-browser-wrap .pb-empty-icon { font-size: 44px; margin-bottom: 14px; }
    .pub-browser-wrap .pb-empty h3 { font-size: 20px; color: #7a90a5; margin-bottom: 6px; }
    .pub-browser-wrap .pb-empty p  { font-size: 13px; color: #9aaabb; }

    /* ---- Responsive ---- */
    @media (max-width: 720px) {
      .pub-browser-wrap .pb-header  { padding: 24px 18px 20px; }
      .pub-browser-wrap .pb-title   { font-size: 24px; }
      .pub-browser-wrap .pb-controls { padding: 10px 14px; }
      .pub-browser-wrap .pb-grid    { padding: 14px; columns: 1; }
    }
  "))
}

# =============================================================================
# MODULE UI
# =============================================================================
#' @param id  A unique string ID for this module instance (e.g. "pubs")
pubBrowserUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$head(
      tags$link(
        href = "https://fonts.googleapis.com/css2?family=Source+Sans+3:wght@300;400;600&display=swap",
        rel  = "stylesheet"
      ),
      .pubBrowserCSS()
    ),
    div(class = "pub-browser-wrap",

      # Header
      div(class = "pb-header",
        div(class = "pb-eyebrow", "Live PubMed Library"),
        div(class = "pb-title",
          tags$em("Anopheles stephensi"), " Research"
        ),
        p(class = "pb-subtitle",
          "Automatically fetched from PubMed · click any card to open the paper")
      ),

      # Controls
      div(class = "pb-controls",
        tags$span(class = "pb-search-icon", "🔍"),
        textInput(ns("search"), label = NULL,
                  placeholder = "Search title, author, journal, keyword, abstract…"),

        div(class = "pb-ctrl-group",
          span(class = "pb-ctrl-label", "From"),
          selectInput(ns("year_from"), label = NULL,
                      choices = c("Any" = "", as.character(2000:2025)),
                      selected = "", width = "85px"),
          span(class = "pb-ctrl-label", "To"),
          selectInput(ns("year_to"), label = NULL,
                      choices = c("Any" = "", as.character(2000:2025)),
                      selected = "", width = "85px")
        ),

        div(class = "pb-ctrl-group",
          span(class = "pb-ctrl-label", "Sort"),
          selectInput(ns("sort"), label = NULL,
                      choices  = c("Newest first" = "desc",
                                   "Oldest first" = "asc",
                                   "A – Z title"  = "alpha"),
                      selected = "desc", width = "128px")
        ),

        actionButton(ns("refresh"), "↻ Refresh"),
        div(class = "pb-count", uiOutput(ns("count")))
      ),

      # Content
      uiOutput(ns("content"))
    )
  )
}

# =============================================================================
# MODULE SERVER
# =============================================================================
#' @param id  Must match the id passed to pubBrowserUI()
pubBrowserServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    pub_data <- reactiveVal(NULL)
    status   <- reactiveVal("loading")

    # force = TRUE bypasses the process-level cache (used by the Refresh button).
    do_fetch <- function(force = FALSE) {
      status("loading")
      pub_data(NULL)
      tryCatch({
        df <- .get_pubmed_papers(force = force)
        if (is.null(df) || nrow(df) == 0) { status("error"); return() }
        pub_data(df)
        status("done")
      }, error = function(e) status("error"))
    }

    # On session start, reuse the cached result if it's still fresh.
    observe({ do_fetch(force = FALSE) })
    observeEvent(input$refresh, { do_fetch(force = TRUE) })

    # Filter + sort
    filtered <- reactive({
      req(status() == "done")
      df <- pub_data()
      if (is.null(df) || nrow(df) == 0) return(df)

      q <- tolower(trimws(input$search))
      if (nchar(q) > 0) {
        df <- df %>% filter(
          str_detect(tolower(title),    fixed(q)) |
          str_detect(tolower(authors),  fixed(q)) |
          str_detect(tolower(keywords), fixed(q)) |
          str_detect(tolower(journal),  fixed(q)) |
          str_detect(tolower(abstract), fixed(q)) |
          str_detect(year, fixed(q))
        )
      }

      yf <- input$year_from
      if (!is.null(yf) && yf != "")
        df <- df %>% filter(!is.na(as.numeric(year)), as.numeric(year) >= as.numeric(yf))

      yt <- input$year_to
      if (!is.null(yt) && yt != "")
        df <- df %>% filter(!is.na(as.numeric(year)), as.numeric(year) <= as.numeric(yt))

      switch(input$sort,
        "desc"  = df %>% arrange(desc(as.numeric(year))),
        "asc"   = df %>% arrange(as.numeric(year)),
        "alpha" = df %>% arrange(title),
        df
      )
    })

    # Count badge
    output$count <- renderUI({
      if (status() != "done") return(NULL)
      n     <- nrow(filtered())
      total <- nrow(pub_data())
      any_filter <- nchar(trimws(input$search)) > 0 ||
                    nchar(input$year_from) > 0 ||
                    nchar(input$year_to)   > 0
      if (any_filter) HTML(paste0("<span>", n, "</span> of ", total, " papers"))
      else            HTML(paste0("<span>", total, "</span> papers"))
    })

    # Main content
    output$content <- renderUI({
      s <- status()

      if (s == "loading") return(
        div(class = "pb-loading",
          div(class = "pb-spinner"),
          p("Fetching papers from PubMed…"),
          tags$small("This usually takes 20 – 40 seconds")
        )
      )

      if (s == "error") return(
        div(class = "pb-error",
          h3("Could not reach PubMed"),
          p("Check your internet connection and click ↻ Refresh to try again.")
        )
      )

      df <- filtered()
      if (is.null(df) || nrow(df) == 0) return(
        div(class = "pb-grid",
          div(class = "pb-empty",
            div(class = "pb-empty-icon", "🦟"),
            h3("No papers found"),
            p("Try adjusting your search or year filters")
          )
        )
      )

      cards <- lapply(seq_len(nrow(df)), function(i) {
        row <- df[i, ]
        kw_list <- if (!is.na(row$keywords) && nchar(row$keywords) > 0) {
          kws <- str_split(row$keywords, ",\\s*")[[1]]
          kws <- kws[seq_len(min(5, length(kws)))]
          lapply(kws, function(k) tags$span(class = "pb-tag", str_trim(k)))
        } else list()

        div(class = "pb-card-wrap",
          tags$a(class = "pb-card", href = row$doi, target = "_blank", rel = "noopener noreferrer",
            div(class = "pb-card-meta",
              span(class = "pb-year",  row$year),
              span(class = "pb-arrow", "↗")
            ),
            div(class = "pb-card-title", row$title),
            div(class = "pb-authors",    row$authors),
            div(class = "pb-journal",    row$journal),
            div(class = "pb-abstract",   row$abstract),
            if (length(kw_list) > 0) div(class = "pb-tags", kw_list)
          )
        )
      })

      div(class = "pb-grid", tagList(cards))
    })
  })
}
