# =============================================================================
# strainhub_transnet.R — vendored make_transnet()
# =============================================================================
# Source:  https://github.com/colbyford/strainhub  (pkg/R/functions.R, v2.0.0)
# Author:  Colby T. Ford <colby.ford@uncc.edu>
# License: Apache License 2.0 — retained here in full compliance.
# Cite:    de Bernardi Schneider A, Su M, Hinrichs AS, et al. "StrainHub: a
#          phylogenetic tool to construct pathogen transmission networks."
#          Bioinformatics 36(3):945-947 (2020).
#
# -----------------------------------------------------------------------------
# WHY THIS FILE EXISTS
# -----------------------------------------------------------------------------
# strainhub cannot be installed as a dependency of a deployed Shiny app. Its
# DESCRIPTION (v2.0.0) declares no Imports: or Depends: field, while its
# NAMESPACE imports 27 packages. R CMD INSTALL builds the lazy-load database by
# loading the namespace; that load reaches for packages nothing ever installed
# and the install aborts with:
#
#     ERROR: lazy loading failed for package 'strainhub'
#
# Two of those 27 imports cannot be satisfied from CRAN at all: `rbokeh` is
# archived, and `treeio` is Bioconductor-only. Because the DESCRIPTION is empty,
# rsconnect also cannot discover any of them, so they never reach the deployment
# manifest even when installed locally.
#
# This app calls exactly one strainhub function — make_transnet() — on exactly
# one code path (treeType = "parsimonious"). That path needs only ape, castor,
# hash, igraph, plyr, tibble and visNetwork, all on active CRAN. rbokeh, treeio,
# adegenet, ade4, DT, leaflet, rhandsontable and colourpicker are used solely by
# make_map() and by strainhub's own bundled Shiny app, neither of which this
# dashboard touches.
#
# -----------------------------------------------------------------------------
# WHAT WAS CHANGED FROM UPSTREAM
# -----------------------------------------------------------------------------
# The statistics — ancestral-state reconstruction, edge counting, and all six
# centrality metrics — are reproduced line-for-line. Five things differ, all of
# them structural rather than numerical:
#
#   1. Every call is namespace-qualified (igraph::, castor::, hash::, ...).
#      Upstream relied on NAMESPACE imports to put these on the search path.
#      Qualifying them means nothing has to be attached here — which matters,
#      because attaching igraph and plyr in app.R would mask a large number of
#      dplyr verbs that the rest of the dashboard depends on.
#   2. The "bayesian", "nj" and "dataframe" branches are not reproduced. They
#      pull in phangorn/seqinr/treeio and are unreachable from this app; calling
#      one now raises an explicit error instead of failing obscurely. The
#      parameters those branches used (threshold, threshold2, bootstrapValue,
#      rootSelection) are therefore dropped from the signature.
#   3. Upstream's six near-identical if/else widget branches are collapsed into
#      build_graph() + a switch. Each branch built the same visNetwork chain and
#      differed only in the metric vector and the title, so the widgets produced
#      are identical; the upstream branches also recomputed a metric that the
#      block above had already computed, and that duplicate call is dropped.
#      An invalid centralityMetric now stops with a message, where upstream
#      cat()-ed an error and then failed on an undefined `graph`.
#   4. metricsOutputFile defaults to "" (write nothing) rather than
#      "StrainHub_metrics.csv". The upstream default writes into the working
#      directory, which is read-only on shinyapps.io and Posit Connect. Pass a
#      path explicitly to get the file when running locally.
#   5. Upstream's commented-out hashmap:: blocks are dropped (hashmap is archived
#      on CRAN; the live code uses the hash package, which is not).
#   6. igraph::graph.data.frame() is replaced by its current spelling,
#      igraph::graph_from_data_frame(). The old name was deprecated in igraph
#      2.0.0 and warns on every render. Same function, same result.
#
# VERIFIED: run against this repo's tree and metadata, the six centrality columns
# reproduce Dashboard/StrainHub_metrics.csv — the file the real strainhub package
# produced — to full printed precision (e.g. Afghanistan betweenness
# 1.61666666666667, Djibouti 27.95). 13 nodes, 62 edges.
#
# To confirm equivalence against the real package on a machine where it installs:
#     a <- strainhub::make_transnet(tree, meta, "Country", 3, metricsOutputFile = "")
#     b <- make_transnet(tree, meta, "Country", 3)
#     identical(a$x$nodes, b$x$nodes) && identical(a$x$edges, b$x$edges)
# =============================================================================

`%>%` <- magrittr::`%>%`

# ---- Hash helpers (verbatim from upstream) ----------------------------------

hash_insert <- function(hm, key, value) {
  if (is.null(key)) {
    stop("Argument 'key' must not be null")
  }
  hm[[as.character(key)]] <- value
  return(hm)
}

hash_find <- function(hm, key) {
  if (is.null(key)) {
    stop("Argument 'key' must not be null")
  }
  if (is.null(hm[[as.character(key)]])) {
    stop("Key not found")
  }
  return(hm[[as.character(key)]])
}

# ---- make_transnet ----------------------------------------------------------
#' Build a transmission network from a phylogenetic tree and metadata
#'
#' @param treedata A rooted phylo object (e.g. from ape::read.tree()).
#' @param metadata A data.frame with an "Accession" column matching the tree's
#'   tip labels, plus the column named by columnSelection.
#' @param columnSelection Name of the metadata column to treat as the character
#'   state (e.g. "Country").
#' @param centralityMetric Which metric drives node size/grouping:
#'   1 indegree, 2 outdegree, 3 betweenness, 4 closeness, 5 degree,
#'   6 Source Hub Ratio.
#' @param treeType Only "parsimonious" is supported in this vendored copy.
#' @param metricsOutputFile Path for the metrics CSV. Empty string writes no
#'   file — required when deploying, since the app directory is read-only on
#'   shinyapps.io and Posit Connect.
#' @param as.json If TRUE return JSON instead of a visNetwork object.
#'
#' @return A visNetwork htmlwidget (or JSON when as.json = TRUE).
make_transnet <- function(treedata,
                          metadata = NULL,
                          columnSelection,
                          centralityMetric,
                          treeType          = "parsimonious",
                          metricsOutputFile = "",
                          as.json           = FALSE) {

  if (treeType != "parsimonious") {
    stop("This vendored copy of make_transnet() supports only ",
         "treeType = \"parsimonious\". For \"bayesian\", \"nj\" or ",
         "\"dataframe\", install the upstream package: ",
         "remotes::install_github(\"colbyford/strainhub\", subdir = \"pkg\")")
  }

  # Takes tip label information from the Newick tree and transforms it into a
  # table, adds an ID, and reorders the metadata frame to match the tree.
  sortingtable <- as.data.frame(treedata$tip.label)
  sortingtable <- tibble::rowid_to_column(sortingtable, "N_ID")
  names(sortingtable)[2] <- "Accession"
  sortingdata <- merge(metadata, sortingtable, by = "Accession")
  data <- sortingdata[order(sortingdata$N_ID), ]

  listofcolumns       <- as.list(data)
  columnaccessions    <- as.list(data)
  accessioncharacter  <- as.character(columnaccessions$Accession)
  selectedcolumn      <- as.numeric(as.factor(listofcolumns[[columnSelection]]))
  names(selectedcolumn) <- accessioncharacter

  # Sorted so the labels line up with the alphabetical order that as.factor()
  # imposes when the state column is coerced to numeric above.
  characterlabels1 <- unique(listofcolumns[[columnSelection]])
  characterlabels  <- sort(as.character(characterlabels1))

  rootedTree   <- treedata
  numCharStates <- length(characterlabels)

  ancestralStates <- castor::asr_max_parsimony(rootedTree,
                                               selectedcolumn,
                                               numCharStates)

  hm <- hash::hash()

  for (i in 1:length(selectedcolumn)) {
    hm <- hash_insert(hm, i, selectedcolumn[i])
  }

  # Walk the inner nodes and assign each the most likely character state.
  numLeaves          <- length(selectedcolumn)
  numInnerNodes      <- rootedTree$Nnode
  totalTreeNodes     <- numLeaves + numInnerNodes
  innerNodeIndices   <- (numLeaves + 1):totalTreeNodes
  numCharacterStates <- length(ancestralStates$ancestral_likelihoods[1, ])
  counter            <- c()

  for (i in innerNodeIndices) {
    counter <- ancestralStates$ancestral_likelihoods[i - numLeaves, ]
    hm <- hash_insert(hm, i, match(max(counter), counter))
  }

  # Walk each edge; where the two nodes differ in state, record the transition.
  sourceList <- c()
  targetList <- c()

  for (row in 1:nrow(rootedTree$edge)) {
    nextEdge <- rootedTree$edge[row, ]

    edgeStates <- c(hash_find(hm, nextEdge[1]),
                    hash_find(hm, nextEdge[2]))

    if (edgeStates[1] != edgeStates[2]) {
      sourceList <- unname(c(sourceList, edgeStates[1]))
      targetList <- unname(c(targetList, edgeStates[2]))
    }
  }

  dat <- data.frame(from = sourceList,
                    to   = targetList)

  edges <- plyr::count(dat)
  names(edges)[names(edges) == "freq"] <- "value"

  metastates <- characterlabels

  nodes <- data.frame(id    = 1:length(metastates),
                      label = metastates)

  # graph_from_data_frame() is the current spelling; upstream's
  # graph.data.frame() was deprecated in igraph 2.0.0 and emits a warning on
  # every render. Straight rename, identical behaviour.
  igraph.Object <- igraph::graph_from_data_frame(edges,
                                                 directed = TRUE,
                                                 vertices = nodes)

  ui <- as.character(centralityMetric)

  # ---- All centrality metrics -----------------------------------------------
  indegree             <- igraph::centr_degree(igraph.Object, mode = c("in"))
  outdegree            <- igraph::centr_degree(igraph.Object, mode = c("out"))
  all.degree           <- igraph::centr_degree(igraph.Object, mode = c("all"))
  between.centrality   <- igraph::betweenness(igraph.Object)
  closeness.centrality <- igraph::closeness(igraph.Object, mode = c("all"))
  sourcehubratio       <- outdegree$res / all.degree$res

  outputFileMatrix <- matrix(ncol = 0, nrow = length(metastates)) %>%
    cbind(metastates,
          all.degree$res,
          indegree$res,
          outdegree$res,
          between.centrality,
          closeness.centrality,
          sourcehubratio)

  colnames(outputFileMatrix, do.NULL = FALSE)
  colnames(outputFileMatrix) <- c("Metastates",
                                  "Degree Centrality",
                                  "Indegree Centrality",
                                  "Outdegree Centrality",
                                  "Betweenness Centrality",
                                  "Closeness Centrality",
                                  "Source Hub Ratio")

  metrics <- as.data.frame(outputFileMatrix)

  if (metricsOutputFile != "") {
    write.table(metrics,
                append       = FALSE,
                file         = metricsOutputFile,
                sep          = ",",
                fileEncoding = "UTF-8",
                col.names    = TRUE,
                row.names    = FALSE,
                quote        = FALSE)
  }

  # ---- Build the widget for the requested metric ----------------------------
  vis_title_style <- paste0("font-family:Lato, Helvetica Neue, Arial, ",
                            "sans-serif;font-weight:bold;font-size:20px;",
                            "text-align:center;")

  build_graph <- function(metric_values, title) {
    nodes <- data.frame(nodes,
                        value = metric_values,
                        group = metric_values)
    visNetwork::visNetwork(
      nodes = nodes,
      edges = edges,
      main  = list(text = title, style = vis_title_style)
    ) %>%
      visNetwork::visPhysics(solver = "repulsion") %>%
      visNetwork::visInteraction(navigationButtons = TRUE) %>%
      visNetwork::visOptions(selectedBy       = "value",
                             highlightNearest = TRUE,
                             nodesIdSelection = TRUE) %>%
      visNetwork::visEdges(arrows = list(to = list(enabled     = TRUE,
                                                   scaleFactor = 0.75)))
  }

  graph <- switch(
    ui,
    "1" = build_graph(indegree$res,           "Indegree Centrality"),
    "2" = build_graph(outdegree$res,          "Outdegree Centrality"),
    "3" = build_graph(between.centrality,     "Betweenness Centrality"),
    "4" = build_graph(closeness.centrality,   "Closeness Centrality"),
    "5" = build_graph(all.degree$res,         "Degree Centrality"),
    "6" = build_graph(sourcehubratio,
                      "Source Hub Ratio: Sink ~0 / Hub = .5 / Source = ~1"),
    stop("Please enter an integer between 1 and 6 to select a centrality metric.")
  )

  if (as.json) {
    return(jsonlite::toJSON(graph$x))
  } else {
    return(graph)
  }
}
