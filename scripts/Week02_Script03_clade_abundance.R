# Week 2 - Task 3
# Match TEsorter Copia/Gypsy clade classifications to the final EDTA
# TE annotation and estimate genomic abundance of each clade.

library(tidyverse)
library(data.table)

# -------------------------------------------------
# Paths
# -------------------------------------------------

WORKDIR <- Sys.getenv(
  "WORKDIR",
  unset = "/data/users/dwalker/organization_annotation_eukaryote_genomes"
)

GFF <- file.path(
  WORKDIR,
  "output/EDTA/assembly.fasta.mod.EDTA.TEanno.gff3"
)

COPIA_CLS <- file.path(
  WORKDIR,
  "output/TEsorter_library/Copia_sequences.fa.rexdb-plant.cls.tsv"
)

GYPSY_CLS <- file.path(
  WORKDIR,
  "output/TEsorter_library/Gypsy_sequences.fa.rexdb-plant.cls.tsv"
)

RESULT_DIR <- file.path(WORKDIR, "results/TE_annotation")
PLOT_DIR <- file.path(WORKDIR, "plots")

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------
# Read TEsorter classification files
# -------------------------------------------------

message("Reading Copia classification...")
copia_cls <- fread(COPIA_CLS, sep = "\t", header = TRUE)
setnames(copia_cls, 1, "TE")

message("Reading Gypsy classification...")
gypsy_cls <- fread(GYPSY_CLS, sep = "\t", header = TRUE)
setnames(gypsy_cls, 1, "TE")

# Remove the EDTA classification suffix:
# TE_00000014#LTR/Copia -> TE_00000014
copia_cls <- copia_cls %>%
  separate(
    TE,
    into = c("Name", "EDTA_classification"),
    sep = "#",
    fill = "right"
  )

gypsy_cls <- gypsy_cls %>%
  separate(
    TE,
    into = c("Name", "EDTA_classification"),
    sep = "#",
    fill = "right"
  )

cls <- bind_rows(copia_cls, gypsy_cls) %>%
  select(
    Name,
    Order,
    Superfamily,
    Clade,
    Complete,
    Strand,
    Domains
  ) %>%
  distinct(Name, .keep_all = TRUE)

message("TEsorter classifications:")
print(table(cls$Superfamily))

# -------------------------------------------------
# Read final EDTA annotation
# -------------------------------------------------

message("Reading EDTA GFF3...")

gff <- read.delim(
  GFF,
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

# Extract EDTA annotation ID and library Name
gff <- gff %>%
  mutate(
    ID = str_match(attributes, "ID=([^;]+)")[, 2],
    Name = str_match(attributes, "Name=([^;]+)")[, 2],
    length = end - start + 1
  )

# Keep Copia and Gypsy annotations only
ltr_gff <- gff %>%
  filter(
    feature %in% c(
      "Copia_LTR_retrotransposon",
      "Gypsy_LTR_retrotransposon"
    )
  )

message("Copia/Gypsy annotation records in final GFF:")
print(table(ltr_gff$feature))

# -------------------------------------------------
# Match EDTA annotations to TEsorter clades
# -------------------------------------------------

matched <- ltr_gff %>%
  inner_join(cls, by = "Name")

message("Matched annotation records:")
print(table(matched$Superfamily))

message("Number of unique matched EDTA library families:")
print(
  matched %>%
    distinct(Name, Superfamily) %>%
    count(Superfamily)
)

# -------------------------------------------------
# Estimate genomic abundance per clade
# -------------------------------------------------

clade_summary <- matched %>%
  group_by(Superfamily, Clade) %>%
  summarise(
    annotation_count = n_distinct(ID),
    library_families = n_distinct(Name),
    bp_annotated = sum(length),
    .groups = "drop"
  ) %>%
  arrange(
    Superfamily,
    desc(annotation_count)
  )

message("")
message("Clade abundance summary:")
print(clade_summary, n = Inf)

# Save full summary
write.table(
  clade_summary,
  file = file.path(
    RESULT_DIR,
    "Week02_Copia_Gypsy_clade_abundance.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------------------------------------
# Save matched annotations for reproducibility
# -------------------------------------------------

matched_output <- matched %>%
  select(
    scaffold,
    start,
    end,
    ID,
    Name,
    Superfamily,
    Clade,
    Complete,
    length
  )

write.table(
  matched_output,
  file = file.path(
    RESULT_DIR,
    "Week02_Copia_Gypsy_matched_annotations.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------------------------------------
# Plot clade abundance
# -------------------------------------------------

plot_data <- clade_summary %>%
  mutate(
    Clade = reorder(Clade, annotation_count)
  )

p <- ggplot(
  plot_data,
  aes(
    x = Clade,
    y = annotation_count,
    fill = Superfamily
  )
) +
  geom_col() +
  coord_flip() +
  facet_wrap(
    ~Superfamily,
    scales = "free_y"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Genomic abundance of Copia and Gypsy clades",
    x = "TEsorter clade",
    y = "Number of EDTA annotation records"
  )

ggsave(
  file.path(PLOT_DIR, "04_Copia_Gypsy_clade_abundance.png"),
  p,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(PLOT_DIR, "04_Copia_Gypsy_clade_abundance.pdf"),
  p,
  width = 10,
  height = 7
)

message("")
message("Finished.")
message("Results: ", RESULT_DIR)
message("Plots:   ", PLOT_DIR)