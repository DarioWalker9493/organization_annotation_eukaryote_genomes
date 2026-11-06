# Week 2 - Task 4
# TE divergence landscape and approximate insertion-age estimation

library(tidyverse)
library(data.table)

# -------------------------------------------------
# Paths
# -------------------------------------------------

WORKDIR <- Sys.getenv(
  "WORKDIR",
  unset = "/data/users/dwalker/organization_annotation_eukaryote_genomes"
)

RNAME <- file.path(
  WORKDIR,
  "output/EDTA/assembly.fasta.mod.EDTA.anno",
  "assembly.fasta.mod.nonzero.out.landscape.Div.Rname.tab"
)

PLOT_DIR <- file.path(WORKDIR, "plots")
RESULT_DIR <- file.path(WORKDIR, "results/TE_annotation")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------
# Substitution rate
# -------------------------------------------------

# Brassicaceae synonymous substitution rate
r <- 8.22e-9

# Approximate age represented by 1% divergence, in million years
age_per_percent_mya <- (0.01 / (2 * r)) / 1e6

message(
  "Approximate age per 1% divergence: ",
  round(age_per_percent_mya, 3),
  " Ma"
)

# -------------------------------------------------
# Read parseRM output
# -------------------------------------------------

message("Reading:")
message(RNAME)

# First two lines contain descriptive/header information.
# Remaining rows contain:
# Rname, Rclass, Rfam and 50 divergence bins.

rep_table <- fread(
  RNAME,
  skip = 2,
  header = FALSE,
  fill = TRUE
)

expected_columns <- 53

if (ncol(rep_table) != expected_columns) {
  stop(
    "Unexpected number of columns in Rname table: ",
    ncol(rep_table),
    " (expected 53)"
  )
}

colnames(rep_table) <- c(
  "Rname",
  "Rclass",
  "Rfam",
  paste0("bin_", 0:49)
)

message("Repeat families read: ", nrow(rep_table))

# -------------------------------------------------
# Convert to long format
# -------------------------------------------------

rep_long <- rep_table %>%
  pivot_longer(
    cols = starts_with("bin_"),
    names_to = "bin",
    values_to = "bp"
  ) %>%
  mutate(
    bin_lower = as.numeric(str_remove(bin, "bin_")),

    # midpoint of each 1%-wide divergence bin
    divergence_percent = bin_lower + 0.5,

    divergence = divergence_percent / 100,

    # T = K / (2r)
    age_years = divergence / (2 * r),
    age_mya = age_years / 1e6,

    fam = paste(Rclass, Rfam, sep = "/")
  )

# -------------------------------------------------
# Follow the course landscape categories
# -------------------------------------------------

rep_long <- rep_long %>%
  filter(
    Rfam != "unknown",
    Rclass != "unknown",
    Rclass %in% c("LTR", "DNA", "MITE"),

    # Course note: exclude Helitrons from this landscape
    fam != "DNA/Helitron"
  )

message("Families retained:")
print(sort(unique(rep_long$fam)))

# -------------------------------------------------
# Sum bp across individual library representatives
# -------------------------------------------------

landscape <- rep_long %>%
  group_by(
    fam,
    Rclass,
    Rfam,
    divergence_percent,
    age_mya
  ) %>%
  summarise(
    bp = sum(bp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Mbp = bp / 1e6
  )

# -------------------------------------------------
# Order families
# -------------------------------------------------

family_order <- landscape %>%
  distinct(fam, Rclass, Rfam) %>%
  mutate(
    group_order = case_when(
      fam == "LTR/Copia" ~ 1,
      fam == "LTR/Gypsy" ~ 2,
      Rclass == "DNA" ~ 3,
      Rclass == "MITE" ~ 4,
      TRUE ~ 5
    )
  ) %>%
  arrange(group_order, Rfam) %>%
  pull(fam)

landscape$fam <- factor(
  landscape$fam,
  levels = unique(family_order)
)

# -------------------------------------------------
# Peak-divergence summary
# -------------------------------------------------

peak_summary <- landscape %>%
  group_by(fam) %>%
  slice_max(
    order_by = bp,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(
    fam,
    peak_divergence_percent = divergence_percent,
    approximate_peak_age_mya = age_mya,
    peak_bp = bp
  ) %>%
  left_join(
    landscape %>%
      group_by(fam) %>%
      summarise(
        total_bp = sum(bp),
        total_Mbp = sum(Mbp),
        .groups = "drop"
      ),
    by = "fam"
  ) %>%
  arrange(desc(total_bp))

message("")
message("Peak divergence / age summary:")
print(peak_summary, n = Inf)

write.table(
  peak_summary,
  file = file.path(
    RESULT_DIR,
    "Week02_TE_landscape_peak_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  landscape,
  file = file.path(
    RESULT_DIR,
    "Week02_TE_divergence_landscape.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------------------------------------
# Plot
# -------------------------------------------------

p <- ggplot(
  landscape,
  aes(
    x = divergence_percent,
    y = Mbp,
    fill = fam
  )
) +
  geom_col(width = 0.95) +
  scale_x_continuous(
    name = "Sequence divergence (%)",
    breaks = seq(0, 50, 5),
    limits = c(0, 50),
    expand = c(0, 0),
    sec.axis = sec_axis(
      ~ . * age_per_percent_mya,
      name = "Approximate insertion age (Ma)"
    )
  ) +
  scale_y_continuous(
    name = "Sequence (Mbp)",
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_brewer(
    palette = "Paired",
    name = "TE superfamily"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  ) +
  labs(
    title = "Transposable element divergence landscape"
  )

ggsave(
  file.path(
    PLOT_DIR,
    "05_TE_divergence_landscape.png"
  ),
  p,
  width = 12,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(
    PLOT_DIR,
    "05_TE_divergence_landscape.pdf"
  ),
  p,
  width = 12,
  height = 7,
  useDingbats = FALSE
)

message("")
message("Finished TE landscape analysis.")
message("Plot directory:   ", PLOT_DIR)
message("Results directory:", RESULT_DIR)