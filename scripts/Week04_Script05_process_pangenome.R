library(tidyverse)


# ============================================================================
# Settings for your GENESPACE run
# ============================================================================

wd <- "/data/users/dwalker/organization_annotation_eukaryote_genomes/output/comparisons/GENESPACE"
focal_genome <- "Lu1"

options(tibble.width = Inf)

# Keep all results from this script together inside the GENESPACE directory.
output_dir <- file.path(wd, "process_pangenome")
plot_dir <- file.path(output_dir, "plots")
table_dir <- file.path(output_dir, "tables")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================================
# STEP 1: Read the TAIR10-reference pan-gene table
# ============================================================================

# pgID = positional pan-gene entry
# genome/id = accession and gene identifier
#
# flag:
#   PASS    = gene assigned through the main orthogroup + synteny grouping
#   array   = additional member of a tandem gene array
#   NSOrtho = non-syntenic orthologue
#
# For the strict syntenic comparison below, only flag == "PASS" is retained.

pg <- read_tsv(
  file.path(wd, "pangenes/TAIR10_pangenes.txt.gz"),
  show_col_types = FALSE
)

message("\nFirst six rows of the pangenome table:")
print(head(pg))

genomes <- sort(unique(pg$genome))
n_genomes <- length(genomes)

message("\nGenomes included in the pangenome table:")
print(genomes)


# Keep unique gene associations for each pan-gene entry.
genes <- pg %>%
  filter(!is.na(id), nzchar(id)) %>%
  distinct(pgID, genome, id)


# ============================================================================
# STEP 2: Count all recorded gene copies per entry and accession
# ============================================================================

counts <- genes %>%
  count(pgID, genome, name = "copies")


# Generate every possible pgID/accession combination so that missing matches
# are represented explicitly as zero copies.
all_combinations <- expand_grid(
  pgID = unique(pg$pgID),
  genome = genomes
)


counts <- all_combinations %>%
  left_join(
    counts,
    by = c("pgID", "genome")
  ) %>%
  mutate(
    copies = replace_na(copies, 0L)
  )


message("\nFirst six pan-gene entry/accession copy counts:")
print(head(counts))


# ============================================================================
# STEP 3: Presence-frequency analysis
# ============================================================================

# Core:
#   entry present in every accession
#
# Accession-specific:
#   entry present in exactly one accession
#
# Accessory:
#   entry present in more than one but fewer than all accessions

entries <- counts %>%
  group_by(pgID) %>%
  summarise(
    n_present = sum(copies > 0),
    single_copy_core = all(copies == 1),
    present_in_focal = any(
      genome == focal_genome & copies > 0
    ),
    .groups = "drop"
  ) %>%
  mutate(
    category = case_when(
      n_present == 0 ~ "no_recorded_match",
      n_present == n_genomes ~ "core",
      n_present == 1 ~ "accession_specific",
      TRUE ~ "accessory"
    )
  )


presence_summary <- entries %>%
  count(category, name = "entries")

message("\nNumber of pan-gene entries in each presence category:")
print(presence_summary)


# A gene can be associated with several positional pan-gene entries.
# Count distinct gene IDs to obtain the number of unique genes represented
# for each accession.

genome_summary <- genes %>%
  distinct(genome, id) %>%
  count(genome, name = "unique_genes")


genome_summary <- expand_grid(
  genome = genomes
) %>%
  left_join(
    genome_summary,
    by = "genome"
  ) %>%
  mutate(
    unique_genes = replace_na(unique_genes, 0L)
  )


message("\nUnique gene IDs represented in the pangenome table:")
print(genome_summary)


# Presence-frequency plot.

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
  geom_col(fill = "steelblue") +
  scale_x_continuous(
    breaks = 0:n_genomes
  ) +
  labs(
    x = "Number of accessions",
    y = "Pan-gene entries",
    title = "Pangenome presence frequency",
    subtitle = "All recorded gene assignments"
  ) +
  theme_minimal()


message("\nPangenome presence-frequency plot:")
print(frequency_plot)


ggsave(
  filename = file.path(
    plot_dir,
    "pangenome_presence_frequency.pdf"
  ),
  plot = frequency_plot,
  width = 7,
  height = 5
)


# ============================================================================
# STEP 4: Count strict syntenic copies using PASS only
# ============================================================================

# The course exercise asks for a comparison between:
#
#   1. All recorded genes
#   2. Syntenic genes only: flag == "PASS"
#
# array and NSOrtho assignments are therefore excluded here.

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
    by = c("pgID", "genome")
  ) %>%
  mutate(
    copies = replace_na(copies, 0L)
  )


# ============================================================================
# STEP 5: Identify copy-variable core entries
# ============================================================================

comparison_counts <- bind_rows(

  counts %>%
    mutate(
      gene_set = "All recorded genes"
    ),

  syntenic_counts %>%
    mutate(
      gene_set = "Syntenic genes (PASS only)"
    )

)


# A copy-variable core entry:
#
#   - is represented in every accession
#   - has different copy numbers among accessions

copy_variable_core <- comparison_counts %>%
  group_by(
    gene_set,
    pgID
  ) %>%
  summarise(
    core = all(copies > 0),
    min_copies = min(copies),
    max_copies = max(copies),
    .groups = "drop"
  ) %>%
  filter(
    core,
    max_copies > min_copies
  ) %>%
  mutate(
    copy_difference = max_copies - min_copies
  )


# Report number of copy-variable core entries in each gene set.
# Include a zero if no entries qualify.

core_summary <- tibble(
  gene_set = unique(comparison_counts$gene_set)
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
    entries = replace_na(entries, 0L)
  )


message(
  paste0(
    "\nNumber of copy-variable core entries in each gene set ",
    "(all recorded genes versus PASS-only syntenic genes):"
  )
)

print(core_summary)


# ============================================================================
# STEP 6: Plot copy-number differences
# ============================================================================

# Select up to 20 copy-variable core entries for each gene set.
#
# Entries are ranked by the largest difference between maximum and minimum
# copy number across the five accessions.

plot_entries <- copy_variable_core %>%
  group_by(gene_set) %>%
  arrange(
    desc(copy_difference),
    desc(max_copies),
    pgID,
    .by_group = TRUE
  ) %>%
  slice_head(n = 20) %>%
  ungroup()


# Create separate row labels so that each facet can contain a different set
# and order of pgIDs.

plot_order <- paste(
  plot_entries$gene_set,
  plot_entries$pgID,
  sep = "::"
)


plot_counts <- comparison_counts %>%
  semi_join(
    plot_entries,
    by = c("gene_set", "pgID")
  ) %>%
  mutate(

    plot_row = factor(
      paste(
        gene_set,
        pgID,
        sep = "::"
      ),
      levels = rev(plot_order)
    ),

    gene_set = factor(
      gene_set,
      levels = core_summary$gene_set
    )

  )


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
        label = copies,
        colour =
          copies >
          min(copies) +
          0.6 * (max(copies) - min(copies))
      ),
      size = 4,
      fontface = "bold"
    ) +

    scale_colour_manual(
      values = c(
        "FALSE" = "white",
        "TRUE" = "black"
      ),
      guide = "none"
    ) +

    facet_wrap(
      ~gene_set,
      scales = "free_y",
      drop = FALSE
    ) +

    scale_y_discrete(
      labels = function(x) {
        sub("^.*::", "", x)
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
    paste0(
      "\nCopy-variable core plot: top entries ranked by ",
      "largest copy-number difference."
    )
  )

  print(cnv_plot)


  ggsave(
    filename = file.path(
      plot_dir,
      "copy_variable_core_comparison.pdf"
    ),
    plot = cnv_plot,
    width = 11,
    height = 9
  )

} else {

  message(
    "No copy-variable core entries were found in either gene set."
  )

}


# ============================================================================
# Save result tables
# ============================================================================

# Presence categories.
write_tsv(
  presence_summary,
  file.path(
    table_dir,
    "pangenome_presence_summary.tsv"
  )
)


# Unique genes represented per accession.
write_tsv(
  genome_summary,
  file.path(
    table_dir,
    "genome_gene_summary.tsv"
  )
)


# Copy-variable core counts.
write_tsv(
  core_summary,
  file.path(
    table_dir,
    "copy_variable_core_summary.tsv"
  )
)


# Detailed copy-variable entries.
write_tsv(
  copy_variable_core,
  file.path(
    table_dir,
    "copy_variable_core_entries.tsv"
  )
)


# List genes associated with pan-gene entries for Lu1 and TAIR10,
# including their GENESPACE assignment flags.

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


message("\n============================================")
message("Pangenome analysis completed successfully")
message("============================================")

message("\nOutput directory:")
message(output_dir)

message("\nPlots:")
message(
  file.path(
    plot_dir,
    "pangenome_presence_frequency.pdf"
  )
)

message(
  file.path(
    plot_dir,
    "copy_variable_core_comparison.pdf"
  )
)

message("\nTables:")
message(table_dir)

