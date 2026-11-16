library(tidyverse)


# ============================================================================
# Settings
# ============================================================================

project_dir <- "/data/users/dwalker/organization_annotation_eukaryote_genomes"

wd <- file.path(
  project_dir,
  "output/comparisons/GENESPACE"
)

focal_genome <- "Lu1"

options(tibble.width = Inf)


# ============================================================================
# Output directories
# ============================================================================

# Analysis tables remain with the GENESPACE results.
output_dir <- file.path(
  wd,
  "process_pangenome"
)

table_dir <- file.path(
  output_dir,
  "tables"
)

# Final presentation-ready plots are stored centrally in the repository.
plot_dir <- file.path(
  project_dir,
  "plots"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================================
# STEP 1: Read the TAIR10-reference pan-gene table
# ============================================================================

# pgID:
#   positional pan-gene entry
#
# genome / id:
#   accession and gene identifier
#
# flag:
#   PASS    = gene assigned through the main orthogroup + synteny grouping
#   array   = additional member of a tandem gene array
#   NSOrtho = non-syntenic orthologue
#
# For the strict syntenic comparison below, only flag == "PASS" is retained.

pangenome_file <- file.path(
  wd,
  "pangenes/TAIR10_pangenes.txt.gz"
)

if (!file.exists(pangenome_file)) {
  stop(
    paste(
      "Missing GENESPACE pangenome file:",
      pangenome_file
    )
  )
}

pg <- read_tsv(
  pangenome_file,
  show_col_types = FALSE
)


message("\n============================================")
message("Week 4 - Pangenome analysis")
message("============================================")


message("\nFirst six rows of the pangenome table:")
print(head(pg))


genomes <- sort(
  unique(pg$genome)
)

n_genomes <- length(genomes)


message("\nGenomes included in the pangenome table:")
print(genomes)


# Keep unique gene associations for each pan-gene entry.

genes <- pg %>%
  filter(
    !is.na(id),
    nzchar(id)
  ) %>%
  distinct(
    pgID,
    genome,
    id
  )


# ============================================================================
# STEP 2: Count all recorded gene copies per entry and accession
# ============================================================================

counts <- genes %>%
  count(
    pgID,
    genome,
    name = "copies"
  )


# Generate every possible pgID/accession combination.
#
# This is important because accessions without a gene assigned to a particular
# pan-gene entry must explicitly receive copy number zero.

all_combinations <- expand_grid(
  pgID = unique(pg$pgID),
  genome = genomes
)


counts <- all_combinations %>%
  left_join(
    counts,
    by = c(
      "pgID",
      "genome"
    )
  ) %>%
  mutate(
    copies = replace_na(
      copies,
      0L
    )
  )


message("\nFirst six pan-gene entry/accession copy counts:")
print(
  head(counts)
)


# ============================================================================
# STEP 3: Presence-frequency analysis
# ============================================================================

# Categories:
#
# core:
#   present in all accessions
#
# accession_specific:
#   present in exactly one accession
#
# accessory:
#   present in more than one but fewer than all accessions

entries <- counts %>%
  group_by(pgID) %>%
  summarise(

    n_present = sum(
      copies > 0
    ),

    single_copy_core = all(
      copies == 1
    ),

    present_in_focal = any(
      genome == focal_genome &
        copies > 0
    ),

    .groups = "drop"

  ) %>%

  mutate(

    category = case_when(

      n_present == 0 ~
        "no_recorded_match",

      n_present == n_genomes ~
        "core",

      n_present == 1 ~
        "accession_specific",

      TRUE ~
        "accessory"

    )

  )


presence_summary <- entries %>%
  count(
    category,
    name = "entries"
  )


message("\nNumber of pan-gene entries in each presence category:")
print(
  presence_summary
)


# ============================================================================
# Unique genes represented per accession
# ============================================================================

# A gene may be associated with several positional pan-gene entries.
#
# Count distinct gene IDs to obtain the actual number of unique genes
# represented for each accession.

genome_summary <- genes %>%
  distinct(
    genome,
    id
  ) %>%
  count(
    genome,
    name = "unique_genes"
  )


genome_summary <- expand_grid(
  genome = genomes
) %>%

  left_join(
    genome_summary,
    by = "genome"
  ) %>%

  mutate(
    unique_genes = replace_na(
      unique_genes,
      0L
    )
  )


message("\nUnique gene IDs represented in the pangenome table:")
print(
  genome_summary
)


# ============================================================================
# Plot 07: Pangenome presence-frequency
# ============================================================================

frequency <- entries %>%
  count(
    n_present,
    name = "pan_gene_entries"
  )


frequency_plot <- ggplot(
  frequency,
  aes(
    x = n_present,
    y = pan_gene_entries
  )
) +

  geom_col(
    fill = "steelblue"
  ) +

  scale_x_continuous(
    breaks = 0:n_genomes
  ) +

  labs(
    x = "Number of accessions",
    y = "Pan-gene entries",
    title = "Pangenome presence frequency",
    subtitle = "All recorded gene assignments"
  ) +

  theme_minimal(
    base_size = 12
  ) +

  theme(
    axis.text = element_text(
      colour = "black"
    ),
    panel.grid.minor = element_blank()
  )


message("\nPlot 07: Pangenome presence frequency")
print(
  frequency_plot
)


ggsave(
  filename = file.path(
    plot_dir,
    "07_pangenome_presence_frequency.pdf"
  ),
  plot = frequency_plot,
  width = 7,
  height = 5
)


ggsave(
  filename = file.path(
    plot_dir,
    "07_pangenome_presence_frequency.png"
  ),
  plot = frequency_plot,
  width = 7,
  height = 5,
  dpi = 300
)


# ============================================================================
# STEP 4: Count strict syntenic copies using PASS only
# ============================================================================

# Course comparison:
#
# 1. All recorded genes
# 2. Syntenic genes only: flag == "PASS"
#
# array and NSOrtho assignments are therefore excluded from the strict
# syntenic comparison.

syntenic_counts <- pg %>%

  filter(
    flag == "PASS",
    !is.na(id),
    nzchar(id)
  ) %>%

  distinct(
    pgID,
    genome,
    id
  ) %>%

  count(
    pgID,
    genome,
    name = "copies"
  )


syntenic_counts <- all_combinations %>%

  left_join(
    syntenic_counts,
    by = c(
      "pgID",
      "genome"
    )
  ) %>%

  mutate(
    copies = replace_na(
      copies,
      0L
    )
  )


# ============================================================================
# STEP 5: Identify copy-variable core entries
# ============================================================================

comparison_counts <- bind_rows(

  counts %>%
    mutate(
      gene_set =
        "All recorded genes"
    ),

  syntenic_counts %>%
    mutate(
      gene_set =
        "Syntenic genes (PASS only)"
    )

)


# Definition:
#
# A copy-variable core entry
#
# - is represented in every accession
# - but has different copy numbers among accessions

copy_variable_core <- comparison_counts %>%

  group_by(
    gene_set,
    pgID
  ) %>%

  summarise(

    core = all(
      copies > 0
    ),

    min_copies = min(
      copies
    ),

    max_copies = max(
      copies
    ),

    .groups = "drop"

  ) %>%

  filter(
    core,
    max_copies > min_copies
  ) %>%

  mutate(
    copy_difference =
      max_copies - min_copies
  )


# ============================================================================
# Summary of copy-variable core entries
# ============================================================================

core_summary <- tibble(
  gene_set =
    unique(
      comparison_counts$gene_set
    )
) %>%

  left_join(

    count(
      copy_variable_core,
      gene_set,
      name = "entries"
    ),

    by = "gene_set"

  ) %>%

  mutate(
    entries = replace_na(
      entries,
      0L
    )
  )


message(
  paste0(
    "\nNumber of copy-variable core entries ",
    "(all recorded genes versus PASS-only syntenic genes):"
  )
)

print(
  core_summary
)


# ============================================================================
# STEP 6: Select copy-variable entries for plotting
# ============================================================================

# For each comparison set, select up to 20 entries with the largest
# copy-number difference across accessions.
#
# Ties are resolved by maximum copy number and then pgID.

plot_entries <- copy_variable_core %>%

  group_by(
    gene_set
  ) %>%

  arrange(
    desc(copy_difference),
    desc(max_copies),
    pgID,
    .by_group = TRUE
  ) %>%

  slice_head(
    n = 20
  ) %>%

  ungroup()


# Different gene sets may contain different pgIDs.
#
# Construct unique plotting labels using both gene_set and pgID.

plot_order <- paste(
  plot_entries$gene_set,
  plot_entries$pgID,
  sep = "::"
)


plot_counts <- comparison_counts %>%

  semi_join(
    plot_entries,
    by = c(
      "gene_set",
      "pgID"
    )
  ) %>%

  mutate(

    plot_row = factor(

      paste(
        gene_set,
        pgID,
        sep = "::"
      ),

      levels = rev(
        plot_order
      )

    ),

    gene_set = factor(
      gene_set,
      levels = core_summary$gene_set
    )

  )


# ============================================================================
# Plot 08: Copy-variable core comparison
# ============================================================================

if (nrow(plot_counts) > 0) {

  cnv_plot <- ggplot(
    plot_counts,
    aes(
      x = genome,
      y = plot_row,
      fill = copies
    )
  ) +

    geom_tile(
      colour = "white"
    ) +

    geom_text(
      aes(
        label = copies
      ),
      size = 4,
      fontface = "bold",
      colour = "black"
    ) +

    facet_wrap(
      ~gene_set,
      scales = "free_y",
      drop = FALSE
    ) +

    scale_y_discrete(
      labels = function(x) {
        sub(
          "^.*::",
          "",
          x
        )
      }
    ) +

    scale_fill_viridis_c() +

    labs(
      x = "Accession",
      y = "Pan-gene entry (pgID)",
      fill = "Copies",
      title = "Copy-variable core entries",
      subtitle =
        "Top 20 per set ranked by largest copy-number difference"
    ) +

    theme_minimal(
      base_size = 12
    ) +

    theme(

      axis.text = element_text(
        colour = "black"
      ),

      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),

      panel.grid = element_blank()

    )


  message(
    "\nPlot 08: Copy-variable core comparison"
  )

  print(
    cnv_plot
  )


  ggsave(
    filename = file.path(
      plot_dir,
      "08_copy_variable_core_comparison.pdf"
    ),
    plot = cnv_plot,
    width = 11,
    height = 9
  )


  ggsave(
    filename = file.path(
      plot_dir,
      "08_copy_variable_core_comparison.png"
    ),
    plot = cnv_plot,
    width = 11,
    height = 9,
    dpi = 300
  )


} else {

  message(
    "No copy-variable core entries were found."
  )

}


# ============================================================================
# Save result tables
# ============================================================================

write_tsv(
  presence_summary,
  file.path(
    table_dir,
    "pangenome_presence_summary.tsv"
  )
)


write_tsv(
  frequency,
  file.path(
    table_dir,
    "pangenome_presence_frequency.tsv"
  )
)


write_tsv(
  genome_summary,
  file.path(
    table_dir,
    "genome_gene_summary.tsv"
  )
)


write_tsv(
  core_summary,
  file.path(
    table_dir,
    "copy_variable_core_summary.tsv"
  )
)


write_tsv(
  copy_variable_core,
  file.path(
    table_dir,
    "copy_variable_core_entries.tsv"
  )
)


write_tsv(
  counts,
  file.path(
    table_dir,
    "all_recorded_copy_counts.tsv"
  )
)


write_tsv(
  syntenic_counts,
  file.path(
    table_dir,
    "PASS_only_copy_counts.tsv"
  )
)


# ============================================================================
# Save Lu1 and TAIR10 gene positions
# ============================================================================

# Useful for later inspection of suspicious high-copy entries.

gene_positions <- pg %>%

  filter(

    genome %in% c(
      focal_genome,
      "TAIR10"
    ),

    !is.na(id),

    nzchar(id)

  ) %>%

  select(
    pgID,
    genome,
    id,
    flag,
    chr,
    start,
    end
  ) %>%

  distinct() %>%

  arrange(
    pgID,
    genome,
    chr,
    start
  )


write_tsv(
  gene_positions,
  file.path(
    table_dir,
    "focal_TAIR10_genes.tsv"
  )
)


# ============================================================================
# Final summary
# ============================================================================

message("\n============================================")
message("Pangenome analysis completed successfully")
message("============================================")


message("\nPresence categories:")
print(
  presence_summary
)


message("\nUnique genes per accession:")
print(
  genome_summary
)


message("\nCopy-variable core entries:")
print(
  core_summary
)


message("\nFinal plots:")

message(
  file.path(
    plot_dir,
    "07_pangenome_presence_frequency.pdf"
  )
)

message(
  file.path(
    plot_dir,
    "07_pangenome_presence_frequency.png"
  )
)

message(
  file.path(
    plot_dir,
    "08_copy_variable_core_comparison.pdf"
  )
)

message(
  file.path(
    plot_dir,
    "08_copy_variable_core_comparison.png"
  )
)


message("\nResult tables:")
message(
  table_dir
)


message("\n============================================")
message("Week 4 pangenome analysis finished")
message("============================================")

