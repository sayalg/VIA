## ---------------------------------------------------------------------------
## VIA: haplotype diversity (Hd) and sampling-effort sensitivity analysis
##
## 1. Reads the COI sequence metadata (Input_Files/genes/COI/metadata.csv).
## 2. Calculates Hd per country (Nei 1987) with standard error and 95% CI,
##    alongside n per country, using the full dataset.
## 3. Repeats the country-by-year subsampling (default 20 per country-year)
##    many times and, for each iteration, rebuilds the haplotype-sharing
##    matrix C = t(H) %*% H and recalculates Hd. Reports edge recovery,
##    edge weight ranges, and the distribution of Hd across iterations.
## 4. Calculates rarefied Hd (equal n per country) as a sensitivity check.
## 5. Tests whether network connectivity tracks sequence count per country.
## 6. Draws the country haplotype-sharing network with edge stability: edges
##    recovered in fewer than 95% of the subsamples are dashed red. Saved as an
##    interactive HTML page plus a CSV of per-edge recovery.
##
## Edit the CONFIG block if needed, then run the whole script.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

## ----------------------------- CONFIG --------------------------------------
## Run from the repository root:  Rscript Scripts/Hapdiversitycalculation.R
csv_path    <- "Input_Files/genes/COI/metadata.csv"  # the dashboard's COI metadata
col_hap     <- "COI"                          # column holding haplotype ID
col_country <- "Country"                      # column holding country
col_year    <- "Year"                         # column holding collection year
n_per_group <- 20                             # sequences per country-year
n_iter      <- 1000                           # number of subsampling repeats
rarefy_n    <- 10                             # equal n per country for rarefied Hd
min_n_hd    <- 2                              # Hd is undefined below this n
seed        <- 123
out_dir     <- "via_sensitivity_output"
## Network with edge stability (section 6)
stable_cut  <- 0.95                           # edges recovered in at least this share of subsamples are solid
base_seed   <- 123                            # seed for the displayed subsample (the dashboard's network)
iter_seed   <- 2026                           # seed for the repeated subsamples
out_html    <- "An_stephensi_Connectivity_with_Recovery.html"
out_csv     <- "network_edge_recovery.csv"
## ---------------------------------------------------------------------------

set.seed(seed)
dir.create(out_dir, showWarnings = FALSE)

## ---- 1. Read and tidy the data --------------------------------------------
## Treat blank cells and the text strings below as missing
na_strings <- c("", "NA", "N/A", "na", "n/a", "NaN", "-", "unknown", "Unknown")
raw <- read.csv(csv_path, na.strings = na_strings, check.names = FALSE)

missing_cols <- setdiff(c(col_hap, col_country, col_year), names(raw))
if (length(missing_cols) > 0) {
  stop("Column(s) not found in the metadata: ", paste(missing_cols, collapse = ", "),
       "\nAvailable columns: ", paste(names(raw), collapse = ", "))
}
dat <- raw |>
  transmute(
    haplotype = as.character(.data[[col_hap]]),
    country   = as.character(.data[[col_country]]),
    year      = as.character(.data[[col_year]])
  ) |>
  filter(!is.na(haplotype), !is.na(country), haplotype != "", country != "")

cat("Records read:", nrow(raw), "\n")
cat("Records used (haplotype and country present):", nrow(dat), "\n")
cat("Records dropped:", nrow(raw) - nrow(dat), "\n")
cat("Countries:", n_distinct(dat$country), "\n")
cat("Haplotypes:", n_distinct(dat$haplotype), "\n\n")

countries <- sort(unique(dat$country))
k <- length(countries)

## ---- 2. Haplotype diversity (Nei 1987) ------------------------------------
hd_nei <- function(haps, min_n = min_n_hd) {
  haps <- haps[!is.na(haps)]
  n <- length(haps)
  if (n < min_n) {
    return(c(n = n, n_hap = length(unique(haps)), Hd = NA_real_, SE = NA_real_))
  }
  p  <- as.numeric(table(haps)) / n
  s2 <- sum(p^2)
  s3 <- sum(p^3)
  Hd <- (n / (n - 1)) * (1 - s2)
  v  <- (2 / (n * (n - 1))) * (2 * (n - 2) * (s3 - s2^2) + s2 - s2^2)
  c(n = n, n_hap = length(p), Hd = Hd, SE = sqrt(max(v, 0)))
}

hd_full <- dat |>
  group_by(country) |>
  summarise(
    res = list(hd_nei(haplotype)),
    n_years = n_distinct(year),
    .groups = "drop"
  ) |>
  mutate(
    n      = vapply(res, `[[`, numeric(1), "n"),
    n_hap  = vapply(res, `[[`, numeric(1), "n_hap"),
    Hd     = vapply(res, `[[`, numeric(1), "Hd"),
    SE     = vapply(res, `[[`, numeric(1), "SE"),
    CI_low = pmax(Hd - 1.96 * SE, 0),
    CI_high = pmin(Hd + 1.96 * SE, 1)
  ) |>
  select(country, n, n_years, n_hap, Hd, SE, CI_low, CI_high) |>
  arrange(desc(n))

write.csv(hd_full, file.path(out_dir, "Hd_full_dataset_by_country.csv"), row.names = FALSE)
cat("Hd on the full dataset (also saved to CSV):\n")
print(as.data.frame(hd_full), digits = 3, row.names = FALSE)
cat("\nNote: when a country has a single haplotype, Hd = 0 and SE = 0, so the\n",
    "interval is uninformative. Interpret low-n countries with caution.\n\n")

## ---- 3. Repeated subsampling ----------------------------------------------
make_H <- function(df) {
  tab <- table(factor(df$haplotype), factor(df$country, levels = countries))
  (unclass(tab) > 0) * 1          # haplotypes x countries, presence/absence
}

subsample_once <- function(df) {
  df |>
    group_by(country, year) |>
    slice_sample(n = n_per_group) |>   # keeps all rows if group has < n_per_group
    ungroup()
}

weights_arr <- array(NA_real_, dim = c(k, k, n_iter),
                     dimnames = list(countries, countries, NULL))
hd_sub      <- matrix(NA_real_, nrow = n_iter, ncol = k, dimnames = list(NULL, countries))
n_sub       <- matrix(NA_real_, nrow = n_iter, ncol = k, dimnames = list(NULL, countries))

cat("Running", n_iter, "subsampling iterations...\n")
for (i in seq_len(n_iter)) {
  sub <- subsample_once(dat)
  H   <- make_H(sub)
  C   <- t(H) %*% H
  diag(C) <- 0
  weights_arr[, , i] <- C
  
  by_c <- split(sub$haplotype, factor(sub$country, levels = countries))
  res  <- vapply(by_c, hd_nei, numeric(4))
  hd_sub[i, ] <- res["Hd", ]
  n_sub[i, ]  <- res["n", ]
}

## Edge recovery and weights for each country pair
pairs <- which(upper.tri(matrix(0, k, k)), arr.ind = TRUE)
edge_summary <- tibble(
  country_a = countries[pairs[, 1]],
  country_b = countries[pairs[, 2]]
) |>
  mutate(
    recovery      = apply(pairs, 1, function(ix) mean(weights_arr[ix[1], ix[2], ] > 0)),
    weight_median = apply(pairs, 1, function(ix) median(weights_arr[ix[1], ix[2], ])),
    weight_min    = apply(pairs, 1, function(ix) min(weights_arr[ix[1], ix[2], ])),
    weight_max    = apply(pairs, 1, function(ix) max(weights_arr[ix[1], ix[2], ]))
  ) |>
  filter(recovery > 0) |>
  arrange(desc(recovery), desc(weight_median))

write.csv(edge_summary, file.path(out_dir, "edge_recovery_across_subsamples.csv"), row.names = FALSE)

## Hd distribution across subsamples
hd_sub_summary <- tibble(
  country        = countries,
  n_subsample_median = apply(n_sub, 2, median),
  Hd_median      = apply(hd_sub, 2, median, na.rm = TRUE),
  Hd_q025        = apply(hd_sub, 2, quantile, probs = 0.025, na.rm = TRUE),
  Hd_q975        = apply(hd_sub, 2, quantile, probs = 0.975, na.rm = TRUE)
) |>
  arrange(desc(n_subsample_median))

write.csv(hd_sub_summary, file.path(out_dir, "Hd_across_subsamples_by_country.csv"), row.names = FALSE)

## ---- 4. Rarefied Hd (equal n per country) ---------------------------------
by_country_full <- split(dat$haplotype, factor(dat$country, levels = countries))
eligible <- names(by_country_full)[vapply(by_country_full, length, integer(1)) >= rarefy_n]

hd_rare <- matrix(NA_real_, nrow = n_iter, ncol = length(eligible),
                  dimnames = list(NULL, eligible))
for (i in seq_len(n_iter)) {
  for (cn in eligible) {
    s <- sample(by_country_full[[cn]], rarefy_n)
    hd_rare[i, cn] <- hd_nei(s)["Hd"]
  }
}

hd_rare_summary <- tibble(
  country   = eligible,
  rarefy_n  = rarefy_n,
  Hd_median = apply(hd_rare, 2, median),
  Hd_q025   = apply(hd_rare, 2, quantile, probs = 0.025),
  Hd_q975   = apply(hd_rare, 2, quantile, probs = 0.975)
) |>
  arrange(desc(Hd_median))

write.csv(hd_rare_summary, file.path(out_dir, "Hd_rarefied_by_country.csv"), row.names = FALSE)
excluded <- setdiff(countries, eligible)
if (length(excluded) > 0) {
  cat("\nCountries with fewer than", rarefy_n, "sequences, excluded from rarefaction:\n  ",
      paste(excluded, collapse = ", "), "\n")
}

## ---- 5. Does connectivity track sequence count? ---------------------------
## Degree = number of countries each country shares at least one haplotype with.
## Total shared = sum of shared-haplotype counts across all partners.
degree_mean <- apply(weights_arr > 0, c(1, 3), sum) |> rowMeans()
shared_mean <- apply(weights_arr, c(1, 3), sum) |> rowMeans()

conn <- tibble(country = countries,
               n_full = hd_full$n[match(countries, hd_full$country)],
               degree_mean = degree_mean[countries],
               shared_total_mean = shared_mean[countries])

cor_degree <- cor.test(conn$n_full, conn$degree_mean, method = "spearman", exact = FALSE)
cor_shared <- cor.test(conn$n_full, conn$shared_total_mean, method = "spearman", exact = FALSE)

cat("\nSpearman correlation, sequences per country vs mean network degree:  rho =",
    round(cor_degree$estimate, 2), " p =", signif(cor_degree$p.value, 3), "\n")
cat("Spearman correlation, sequences per country vs mean total shared haplotypes: rho =",
    round(cor_shared$estimate, 2), " p =", signif(cor_shared$p.value, 3), "\n")

write.csv(conn, file.path(out_dir, "connectivity_vs_sample_size.csv"), row.names = FALSE)

## ---- 6. Country network with edge stability --------------------------------
## The country haplotype-sharing network, with each edge styled by how often it
## reappears across repeated subsamples: solid grey edges were recovered in at
## least `stable_cut` of them, dashed red edges in fewer. Uses the tidy data
## above and the CONFIG settings. The repeats here use `iter_seed`, so the
## recovery values differ slightly from section 3 (which uses `seed`); the
## network figure uses these.
library(visNetwork)
library(igraph)

shared_matrix <- function(df) {
  H <- make_H(df)
  C <- t(H) %*% H                  # shared unique haplotypes between countries
  diag(C) <- 0
  C
}

## The displayed subsample: same seed and steps as the dashboard's network.
set.seed(base_seed)
dat_displayed <- subsample_once(dat)
base_matrix   <- shared_matrix(dat_displayed)

## Repeated subsamples, with their own seed.
set.seed(iter_seed)
net_weights <- array(0, dim = c(k, k, n_iter),
                     dimnames = list(countries, countries, NULL))
cat("Running", n_iter, "subsampling iterations for edge stability...\n")
for (i in seq_len(n_iter)) {
  net_weights[, , i] <- shared_matrix(subsample_once(dat))
}

net_edges <- data.frame(
  country_a     = countries[pairs[, 1]],
  country_b     = countries[pairs[, 2]],
  base_weight   = base_matrix[pairs],
  recovery      = apply(net_weights, c(1, 2), function(x) mean(x > 0))[pairs],
  weight_median = apply(net_weights, c(1, 2), median)[pairs],
  weight_min    = apply(net_weights, c(1, 2), min)[pairs],
  weight_max    = apply(net_weights, c(1, 2), max)[pairs],
  stringsAsFactors = FALSE
) |>
  filter(base_weight > 0 | recovery > 0) |>
  arrange(desc(recovery), desc(weight_median))

write.csv(net_edges, file.path(out_dir, out_csv), row.names = FALSE)
cat("Edge table written to", file.path(out_dir, out_csv), "\n")

## Nodes: size = sequences in the displayed subsample.
n_full <- dat           |> count(country, name = "n_full")
n_sub  <- dat_displayed |> count(country, name = "n_sub")

net_nodes <- data.frame(country = countries, stringsAsFactors = FALSE) |>
  left_join(n_full, by = "country") |>
  left_join(n_sub,  by = "country") |>
  mutate(id    = country,
         label = paste0(country, " (n = ", n_sub, ")"),
         group = country,
         value = n_sub,
         title = paste0("<b>", country, "</b><br>Sequences in displayed subsample: ", n_sub,
                        "<br>Sequences in full dataset: ", n_full))

## Edges shown are those in the displayed subsample, styled by stability.
net_links <- net_edges |>
  filter(base_weight > 0) |>
  transmute(
    from   = country_a,
    to     = country_b,
    value  = base_weight,                    # thickness = shared haplotypes
    dashes = recovery < stable_cut,
    color  = ifelse(recovery < stable_cut, "rgba(200,60,40,0.85)", "rgba(100,100,100,0.45)"),
    title  = paste0("Shared haplotypes (displayed subsample): ", base_weight,
                    "<br>Median across subsamples: ", weight_median,
                    " (range ", weight_min, " to ", weight_max, ")",
                    "<br>Recovered in ", round(100 * recovery, 1), "% of subsamples")
  )

## Node colours by country (both spellings of Saudi Arabia, as the data and the
## older figures differ).
node_colours <- c(Afghanistan = "#CD9B1D", Djibouti = "#00008B", Ethiopia = "#009ACD",
                  India = "#EE9A00", Iran = "#CD853F", Kenya = "#6CA6CD", Niger = "#7CCD7C",
                  Pakistan = "#EEDC82", SaudiArabia = "#556B2F", "Saudi Arabia" = "#556B2F",
                  Somalia = "#BCD2EE", "Sri Lanka" = "#BFEFFF", Sudan = "#1874CD",
                  "United Arab Emirates" = "#FF7F00", Yemen = "#4169E1")

legend_edges <- data.frame(
  label  = c(paste0("Stable (>= ", 100 * stable_cut, "%)"),
             paste0("Less stable (< ", 100 * stable_cut, "%)")),
  color  = c("rgba(100,100,100,0.6)", "rgba(200,60,40,0.9)"),
  dashes = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)

net_vis <- visNetwork(net_nodes, net_links,
                      main    = "An. stephensi: Country Connectivity",
                      submain = paste0("Edges show shared haplotypes in one random subsample; ",
                                       "dashed red edges appeared in fewer than ",
                                       100 * stable_cut, "% of ", n_iter, " subsamples")) |>
  visIgraphLayout(layout = "layout_with_fr", randomSeed = 123) |>
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) |>
  visEdges(smooth = TRUE) |>
  visPhysics(solver = "forceAtlas2Based", stabilization = TRUE) |>
  visLegend(addEdges = legend_edges, useGroups = FALSE, addNodes = NULL,
            main = "Edge stability", zoom = FALSE)

for (g in names(node_colours)) {
  net_vis <- visGroups(net_vis, groupname = g, color = node_colours[[g]])
}

visSave(net_vis, file = file.path(out_dir, out_html))
cat("Network saved to", file.path(out_dir, out_html), "\n")

cat("\nDone. Outputs written to:", normalizePath(out_dir), "\n")