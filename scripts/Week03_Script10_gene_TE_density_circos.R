# Week 3
# Compare final MAKER gene density with EDTA TE density
# across the 10 largest Lu-1 scaffolds.

library(tidyverse)
library(data.table)
library(circlize)

# -------------------------------------------------
# Paths
# -------------------------------------------------

WORKDIR <- Sys.getenv(
  "WORKDIR",
  unset = "/data/users/dwalker/organization_annotation_eukaryote_genomes"
)

te_gff_file <- file.path(
  WORKDIR,
  "output/EDTA/assembly.fasta.mod.EDTA.TEanno.gff3"
)

gene_gff_file <- file.path(
  WORKDIR,
  "output/MAKER/final/filtered.genes.renamed.gff3"
)

fai_file <- file.path(
  WORKDIR,
  "data/assemblies/Lu1_Flye.fasta.fai"
)

plot_dir <- file.path(WORKDIR, "plots")
result_dir <- file.path(WORKDIR, "results/gene_annotation")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------
# Read scaffold lengths
# -------------------------------------------------

fai <- fread(
  fai_file,
  header = FALSE,
  select = 1:2,
  col.names = c("scaffold", "length")
)

n_scaffolds <- 10

top_scaffolds <- fai %>%
  arrange(desc(length)) %>%
  slice_head(n = n_scaffolds)

message("Largest scaffolds:")
print(top_scaffolds)

ideogram <- top_scaffolds %>%
  transmute(
    scaffold = scaffold,
    start = 0,
    end = length
  )

# -------------------------------------------------
# GFF reader
# -------------------------------------------------

read_gff <- function(path) {

  x <- read.delim(
    path,
    sep = "\t",
    header = FALSE,
    comment.char = "#",
    quote = "",
    stringsAsFactors = FALSE
  )

  colnames(x) <- c(
    "scaffold",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes"
  )

  x
}

# -------------------------------------------------
# Final MAKER genes
# -------------------------------------------------

message("Reading final MAKER gene annotation...")

gene_gff <- read_gff(gene_gff_file)

genes <- gene_gff %>%
  filter(
    feature == "gene",
    scaffold %in% top_scaffolds$scaffold
  ) %>%
  mutate(
    start0 = start - 1,
    midpoint = floor((start0 + end) / 2)
  )

message("Genes on selected scaffolds: ", nrow(genes))

gene_bed <- genes %>%
  select(scaffold, start0, end)

# -------------------------------------------------
# EDTA TEs
# -------------------------------------------------

message("Reading EDTA annotation...")

te_gff <- read_gff(te_gff_file)

te <- te_gff %>%
  mutate(
    TE_group = case_when(
      feature == "Gypsy_LTR_retrotransposon" ~ "Gypsy",
      feature == "Copia_LTR_retrotransposon" ~ "Copia",
      grepl("_TIR_transposon$", feature) ~ "TIR",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(TE_group),
    scaffold %in% top_scaffolds$scaffold
  ) %>%
  mutate(
    start0 = start - 1,
    midpoint = floor((start0 + end) / 2)
  )

gypsy <- te %>%
  filter(TE_group == "Gypsy") %>%
  select(scaffold, start0, end)

copia <- te %>%
  filter(TE_group == "Copia") %>%
  select(scaffold, start0, end)

tir <- te %>%
  filter(TE_group == "TIR") %>%
  select(scaffold, start0, end)

# -------------------------------------------------
# Scaffold-level summary
# -------------------------------------------------

gene_summary <- genes %>%
  count(scaffold, name = "feature_count") %>%
  complete(
    scaffold = top_scaffolds$scaffold,
    fill = list(feature_count = 0)
  ) %>%
  mutate(feature_group = "Gene")

te_summary <- te %>%
  count(scaffold, TE_group, name = "feature_count") %>%
  complete(
    scaffold = top_scaffolds$scaffold,
    TE_group = c("Gypsy", "Copia", "TIR"),
    fill = list(feature_count = 0)
  ) %>%
  rename(feature_group = TE_group)

summary_table <- bind_rows(
  gene_summary,
  te_summary
) %>%
  left_join(top_scaffolds, by = "scaffold") %>%
  mutate(
    features_per_Mb = feature_count / (length / 1e6)
  ) %>%
  arrange(
    factor(scaffold, levels = top_scaffolds$scaffold),
    feature_group
  )

write.table(
  summary_table,
  file = file.path(
    result_dir,
    "Week03_gene_TE_density_top10_scaffolds.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------------------------------------
# 100-kb windows for quantitative comparison
# -------------------------------------------------

window_size <- 100000

windows <- map2_dfr(
  top_scaffolds$scaffold,
  top_scaffolds$length,
  function(scaffold_name, scaffold_length) {

    starts <- seq(
      0,
      scaffold_length - 1,
      by = window_size
    )

    tibble(
      scaffold = scaffold_name,
      window_start = starts,
      window_end = pmin(
        starts + window_size,
        scaffold_length
      )
    )
  }
)

gene_windows <- genes %>%
  mutate(
    window_start =
      floor(midpoint / window_size) * window_size
  ) %>%
  count(
    scaffold,
    window_start,
    name = "Gene"
  )

te_windows <- te %>%
  mutate(
    window_start =
      floor(midpoint / window_size) * window_size
  ) %>%
  count(
    scaffold,
    window_start,
    TE_group
  ) %>%
  pivot_wider(
    names_from = TE_group,
    values_from = n,
    values_fill = list(n = 0)
  )

window_table <- windows %>%
  left_join(
    gene_windows,
    by = c("scaffold", "window_start")
  ) %>%
  left_join(
    te_windows,
    by = c("scaffold", "window_start")
  ) %>%
  mutate(
    across(
      c(Gene, Gypsy, Copia, TIR),
      ~ replace_na(.x, 0)
    ),
    Total_TE = Gypsy + Copia + TIR,
    window_Mb =
      (window_end - window_start) / 1e6,
    Gene_per_Mb = Gene / window_Mb,
    Gypsy_per_Mb = Gypsy / window_Mb,
    Copia_per_Mb = Copia / window_Mb,
    TIR_per_Mb = TIR / window_Mb,
    Total_TE_per_Mb = Total_TE / window_Mb
  )

write.table(
  window_table,
  file = file.path(
    result_dir,
    "Week03_gene_TE_density_100kb_windows.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------------------------------------
# Gene-vs-TE correlations
# -------------------------------------------------

correlation_table <- tibble(
  TE_group = c(
    "Gypsy",
    "Copia",
    "TIR",
    "Total_TE"
  ),
  spearman_rho = c(
    cor(
      window_table$Gene_per_Mb,
      window_table$Gypsy_per_Mb,
      method = "spearman"
    ),
    cor(
      window_table$Gene_per_Mb,
      window_table$Copia_per_Mb,
      method = "spearman"
    ),
    cor(
      window_table$Gene_per_Mb,
      window_table$TIR_per_Mb,
      method = "spearman"
    ),
    cor(
      window_table$Gene_per_Mb,
      window_table$Total_TE_per_Mb,
      method = "spearman"
    )
  )
)

write.table(
  correlation_table,
  file = file.path(
    result_dir,
    "Week03_gene_TE_density_correlations.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Gene-vs-TE Spearman correlations:")
print(correlation_table)

# -------------------------------------------------
# Circos plot
# -------------------------------------------------

make_circos_plot <- function() {

  circos.clear()

  circos.par(
    start.degree = 90,
    gap.degree = 3,
    track.margin = c(0.01, 0.01)
  )

  circos.genomicInitialize(
    ideogram,
    plotType = c("axis", "labels"),
    major.by = 1000000
  )

  # Outer track: genes
  circos.genomicDensity(
    gene_bed,
    window.size = window_size,
    col = "#000000",
    track.height = 0.10
  )

  # Gypsy
  circos.genomicDensity(
    gypsy,
    window.size = window_size,
    col = "#D55E00",
    track.height = 0.10
  )

  # Copia
  circos.genomicDensity(
    copia,
    window.size = window_size,
    col = "#009E73",
    track.height = 0.10
  )

  # Inner track: TIR
  circos.genomicDensity(
    tir,
    window.size = window_size,
    col = "#0072B2",
    track.height = 0.10
  )

  legend(
    "topleft",
    legend = c(
      "Genes",
      "Gypsy",
      "Copia",
      "TIR"
    ),
    fill = c(
      "#000000",
      "#D55E00",
      "#009E73",
      "#0072B2"
    ),
    border = NA,
    bty = "n",
    cex = 0.85
  )

  title(
    "Gene and TE density across the 10 largest Lu-1 scaffolds",
    cex.main = 1.1
  )

  circos.clear()
}

# -------------------------------------------------
# Save PNG
# -------------------------------------------------

png(
  file.path(
    plot_dir,
    "06_Gene_TE_density_circos.png"
  ),
  width = 3000,
  height = 3000,
  res = 300
)

make_circos_plot()
dev.off()

# -------------------------------------------------
# Save PDF
# -------------------------------------------------

pdf(
  file.path(
    plot_dir,
    "06_Gene_TE_density_circos.pdf"
  ),
  width = 10,
  height = 10
)

make_circos_plot()
dev.off()

message("Finished.")
message("Plots: ", plot_dir)
message("Results: ", result_dir)
