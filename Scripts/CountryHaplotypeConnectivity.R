# 1. Load Libraries
library(readxl)
library(visNetwork)
library(tidyverse)
library(igraph)

# 2. Set Workspace and Load Data
# Ensure these filenames match your folder exactly
setwd("C:/Users/Waymi/Box/Carter Lab/Projects/NIH Surveillance")
metadata <- read_excel("FinalCOI_metadata.xlsx")
dna <- read.dna("COI_representative317.fasta", format = "fasta")

#  3. Clean Metadata 
# Remove any rows where Country or Haplotype is missing to prevent ghost nodes
metadata_clean <- metadata %>%
  filter(!is.na(Country) & Country != "" & Country != "NA") %>%
  filter(!is.na(Haplotype) & Haplotype != "" & Haplotype != "NA")

# 4.Random Subsampling
set.seed(123) 
metadata_subsampled <- metadata_clean %>%
  group_by(Country, Year) %>%
  # Randomly pick 20 samples per group. If a group has < 20, it takes all of them.
  slice_sample(n = 20) %>% 
  ungroup()

# 5. Binary Haplotype Logic
# We create a table of Haplotypes vs Country
h_matrix_raw <- table(metadata_subsampled$Haplotype, metadata_subsampled$Country)

# 6. Convert to Binary: 1 if the country has the haplotype, 0 if not.
h_matrix_binary <- h_matrix_raw
h_matrix_binary[h_matrix_binary > 0] <- 1

# 7. Generate the shared unique haplotype count
country_matrix <- as.matrix(t(h_matrix_binary) %*% h_matrix_binary)
diag(country_matrix) <- 0

# 8. Define Nodes (Countries)
country_totals <- metadata_subsampled %>%
  group_by(Country) %>%
  summarize(n = n())

nodes_country <- data.frame(
  id = 1:nrow(country_totals),
  label = country_totals$Country,
  group = country_totals$Country,
  value = country_totals$n,        # Node size = total samples in the 20-cap set
  title = paste("Included Samples:", country_totals$n)
)

# 6. Define Links (Shared Unique Haplotypes) 
edges_idx <- which(country_matrix > 0, arr.ind = TRUE)
edges_idx <- edges_idx[edges_idx[,1] < edges_idx[,2], ]

links_country <- data.frame(
  from = edges_idx[, 1],
  to = edges_idx[, 2],
  # Thickness is the count of shared unique Haplotype Names
  value = mapply(function(f, t) country_matrix[f, t], edges_idx[, 1], edges_idx[, 2])
)
# Note: We are NOT adding a 'label' column here, so no numbers will appear on lines.

# 7. Final Visualization with Phylogenetic Colors
vis_country_fair <- visNetwork(nodes_country, links_country, 
                               main = "An. stephensi: Country Connectivity") %>%
  visIgraphLayout(layout = "layout_with_fr") %>% 
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visEdges(color = "rgba(100, 100, 100, 0.4)", smooth = TRUE) %>%
  # Your specific Color Palette
  visGroups(groupname = "Afghanistan", color = "#CD9B1D") %>%
  visGroups(groupname = "Djibouti", color = "#00008B") %>%
  visGroups(groupname = "Ethiopia", color = "#009ACD") %>%
  visGroups(groupname = "India", color = "#EE9A00") %>%
  visGroups(groupname = "Iran", color = "#CD853F") %>%
  visGroups(groupname = "Kenya", color = "#6CA6CD") %>%
  visGroups(groupname = "Niger", color = "#7CCD7C") %>%
  visGroups(groupname = "Pakistan", color = "#EEDC82") %>%
  visGroups(groupname = "SaudiArabia", color = "#556B2F") %>% 
  visGroups(groupname = "Somalia", color = "#BCD2EE") %>%
  visGroups(groupname = "Sri Lanka", color = "#BFEFFF") %>%
  visGroups(groupname = "Sudan", color = "#1874CD") %>%
  visGroups(groupname = "United Arab Emirates", color = "#FF7F00") %>%
  visGroups(groupname = "Yemen", color = "#4169E1") %>%
  visPhysics(solver = "forceAtlas2Based", stabilization = TRUE)

# Save
visSave(vis_country_fair, file = "An_stephensi_Simplified_Connectivity.html")
