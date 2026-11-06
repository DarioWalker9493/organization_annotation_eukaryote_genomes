# Week 1 - Task 2.2
# Visualize the genomic distribution of Gypsy, Copia and TIR elements
# across the largest scaffolds of the Lu-1 Flye assembly.

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

gff_file <- file.path(
  WORKDIR,
  "output/EDTA/assembly.fasta.mod.EDTA.TEanno.gff3"
)

fai_file <- file.path(
  WORKDIR,
  "data/assemblies/Lu1_Flye.fasta.fai"
)

plot_dir <- file.path(WORKDIR, "plots")
result_dir <- file.path(WORKDIR, "results/TE_annotation")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------
# Read scaffold lengths
# -------------------------------------------------

message("Reading assembly index: ", fai_file)

fai <- fread(
  fai_file,
  header = FALSE,
  select = 1:2,
  col.names = c("scaffold", "length")
)

# Use the 10 largest scaffolds for visualization.
n_scaffolds <- 10

top_scaffolds <- fai %>%
  arrange(desc(length)) %>%
  slice_head(n = n_scaffolds)

message("Largest scaffolds:")
print(top_scaffolds)

# Custom ideogram required by circlize:
# scaffold, start, end
ideogram <- top_scaffolds %>%
  transmute(
    scaffold = scaffold,
    start = 0,
    end = length
  )

# -------------------------------------------------
# Read EDTA annotation
# -------------------------------------------------

message("Reading EDTA annotation: ", gff_file)

gff <- read.delim(
  gff_file,
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  quote = "",
  stringsAsFactors = FALSE
)

colnames(gff) <- c(
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

# -------------------------------------------------
# Classify the TE groups needed for Task 2.2
# -------------------------------------------------

te <- gff %>%
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
    # Convert GFF start coordinate to zero-based start for genomic plotting.
    start0 = start - 1
  )

message("TE counts on selected scaffolds:")
print(table(te$TE_group))

# -------------------------------------------------
# Summary table
# -------------------------------------------------

summary_table <- te %>%
  count(scaffold, TE_group, name = "element_count") %>%
  complete(
    scaffold = top_scaffolds$scaffold,
    TE_group = c("Gypsy", "Copia", "TIR"),
    fill = list(element_count = 0)
  ) %>%
  left_join(top_scaffolds, by = "scaffold") %>%
  mutate(
    elements_per_Mb = element_count / (length / 1e6)
  ) %>%
  arrange(
    factor(scaffold, levels = top_scaffolds$scaffold),
    TE_group
  )

write.table(
  summary_table,
  file = file.path(
    result_dir,
    "Week02_TE_density_top10_scaffolds.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------------------------------------
# Prepare BED-like data for circlize
# -------------------------------------------------

gypsy <- te %>%
  filter(TE_group == "Gypsy") %>%
  select(scaffold, start0, end)

copia <- te %>%
  filter(TE_group == "Copia") %>%
  select(scaffold, start0, end)

tir <- te %>%
  filter(TE_group == "TIR") %>%
  select(scaffold, start0, end)

# 100 kb density windows
window_size <- 100000

# -------------------------------------------------
# Circos plotting function
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

  # Outer track: Gypsy
  circos.genomicDensity(
    gypsy,
    window.size = window_size,
    col = "#D55E00",
    track.height = 0.12
  )

  # Middle track: Copia
  circos.genomicDensity(
    copia,
    window.size = window_size,
    col = "#009E73",
    track.height = 0.12
  )

  # Inner track: TIR
  circos.genomicDensity(
    tir,
    window.size = window_size,
    col = "#0072B2",
    track.height = 0.12
  )

  legend(
    "topleft",
    legend = c("Gypsy", "Copia", "TIR"),
    fill = c("#D55E00", "#009E73", "#0072B2"),
    border = NA,
    bty = "n",
    cex = 0.9
  )

  title(
    "TE density across the 10 largest Lu-1 scaffolds",
    cex.main = 1.2
  )

  circos.clear()
}

# -------------------------------------------------
# Save PNG
# -------------------------------------------------

png(
  file.path(plot_dir, "03_TE_density_circos.png"),
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
  file.path(plot_dir, "03_TE_density_circos.pdf"),
  width = 10,
  height = 10
)

make_circos_plot()

dev.off()

message("Finished.")
message("Plots written to: ", plot_dir)
message("Summary written to: ", result_dir)