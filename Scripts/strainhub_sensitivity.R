# Sampling sensitivity of the StrainHub phylogeographic dispersal network.
#
# Repeats the dashboard's network inference (Tools/build_app_cache.R ->
# make_transnet(): maximum-parsimony ASR of Country with castor, equal costs)
# on balanced random subsamples of tips from the same COI ML tree, then reports
# how often each edge of the full network is recovered and how stable
# betweenness centrality is.
#
# Run from the repository root:  Rscript Scripts/strainhub_sensitivity.R
# ponytail: tips are resampled on the fixed ML tree; re-running RAxML per
# subsample would also capture topological uncertainty.

source("Tools/strainhub_transnet.R")

N_PER_COUNTRY <- 10    # countries with fewer tips keep all of theirs
N_REPS        <- 100
set.seed(1)

tree <- ape::read.tree("Input_Files/genes/COI/tree.newick")
md   <- read.csv("Input_Files/genes/COI/metadata.csv")
md   <- data.frame(Accession = md$FASTA.ID, Country = md$Country)
md   <- md[md$Accession %in% tree$tip.label & !is.na(md$Country) & md$Country != "", ]

net <- function(tips) {
  g <- make_transnet(ape::keep.tip(tree, tips), md[md$Accession %in% tips, ],
                     columnSelection = "Country", centralityMetric = 3,
                     metricsOutputFile = "")$x
  lab <- setNames(g$nodes$label, g$nodes$id)
  list(edges = paste(lab[as.character(g$edges$from)], "->", lab[as.character(g$edges$to)]),
       betw  = setNames(g$nodes$value, g$nodes$label))
}

full <- net(md$Accession)
countries <- names(full$betw)

reps <- replicate(N_REPS, simplify = FALSE, {
  tips <- unlist(lapply(split(md$Accession, md$Country),
                        function(x) if (length(x) > N_PER_COUNTRY) sample(x, N_PER_COUNTRY) else x))
  net(tips)
})

edge_support <- sapply(full$edges, function(e) mean(sapply(reps, function(r) e %in% r$edges)))
new_edges    <- table(unlist(lapply(reps, function(r) setdiff(r$edges, full$edges))))
rho <- sapply(reps, function(r) cor(full$betw, r$betw[countries], method = "spearman"))
top <- names(which.max(full$betw))
top_kept <- mean(sapply(reps, function(r) names(which.max(r$betw)) == top))

cat(sprintf("Full network: %d tips, %d countries, %d directed edges\n",
            nrow(md), length(countries), length(full$edges)))
cat(sprintf("Replicates: %d, up to %d tips per country\n\n", N_REPS, N_PER_COUNTRY))
cat("Recovery of full-network edges (fraction of replicates):\n")
print(round(sort(edge_support, decreasing = TRUE), 2))
cat(sprintf("\nMedian edge recovery: %.2f; edges recovered in >=50%% of replicates: %d/%d\n",
            median(edge_support), sum(edge_support >= 0.5), length(edge_support)))
cat(sprintf("Edges absent from the full network appearing in >=50%% of replicates: %d\n",
            sum(new_edges >= N_REPS / 2)))
cat(sprintf("Betweenness vs full network, Spearman rho: median %.2f (IQR %.2f-%.2f)\n",
            median(rho, na.rm = TRUE), quantile(rho, .25, na.rm = TRUE), quantile(rho, .75, na.rm = TRUE)))
cat(sprintf("Highest-betweenness country in full network (%s) also highest in %.0f%% of replicates\n",
            top, 100 * top_kept))

dir.create("Figures", showWarnings = FALSE)
write.csv(data.frame(edge = names(edge_support), recovery = unname(edge_support)),
          "Figures/strainhub_edge_recovery.csv", row.names = FALSE)
