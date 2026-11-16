library(GENESPACE)


# ============================================================================
# Settings
# ============================================================================

project_dir <- "/data/users/dwalker/organization_annotation_eukaryote_genomes"

wd <- file.path(
  project_dir,
  "output/comparisons/GENESPACE"
)

focal_genome <- "Lu1"

plot_dir <- file.path(
  project_dir,
  "plots"
)

riparian_dir <- file.path(
  wd,
  "processed_riparian"
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  riparian_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================================
# Load completed GENESPACE run
# ============================================================================

gsparams_file <- file.path(
  wd,
  "results",
  "gsParams.rda"
)

if (!file.exists(gsparams_file)) {
  stop(
    paste(
      "Missing GENESPACE parameter file:",
      gsparams_file
    )
  )
}

load(gsparams_file)

if (!exists("gsParam")) {
  stop("Object gsParam was not found in gsParams.rda")
}


message("============================================")
message("Week 5 - GENESPACE synteny visualization")
message("============================================")

message("")
message("GENESPACE directory:")
message(wd)

message("")
message("Focal genome:")
message(focal_genome)


# ============================================================================
# PLOT 09
# Lu1 versus TAIR10
# ============================================================================

pairwise_pdf <- file.path(
  plot_dir,
  "09_Lu1_vs_TAIR10_riparian.pdf"
)

message("")
message("Creating Plot 09:")
message(pairwise_pdf)


plot_riparian(
  gsParam,
  genomeIDs = c(
    "TAIR10",
    focal_genome
  ),
  refGenome = "TAIR10",
  reorderBySynteny = TRUE,
  pdfFile = pairwise_pdf
)


# Also keep an analysis copy with the GENESPACE results.

file.copy(
  pairwise_pdf,
  file.path(
    riparian_dir,
    "Lu1_vs_TAIR10.pdf"
  ),
  overwrite = TRUE
)


# ============================================================================
# PLOT 10
# All five genomes relative to TAIR10
# ============================================================================

all_genomes_pdf <- file.path(
  plot_dir,
  "10_all_genomes_TAIR10_riparian.pdf"
)

message("")
message("Creating Plot 10:")
message(all_genomes_pdf)


plot_riparian(
  gsParam,
  refGenome = "TAIR10",
  reorderBySynteny = TRUE,
  syntenyWeight = 1,
  scalePlotWidth = 1.5,
  pdfFile = all_genomes_pdf
)


file.copy(
  all_genomes_pdf,
  file.path(
    riparian_dir,
    "all_genomes_TAIR10.pdf"
  ),
  overwrite = TRUE
)


# ============================================================================
# Orientation note
# ============================================================================

# The course exercise allows individual Lu1 scaffolds to be flipped for
# DISPLAY if their main syntenic connection runs in reverse.
#
# We deliberately do not define invertTheseChrs yet.
#
# First inspect Plot 09.
# If a Lu1 scaffold clearly shows a reversed main syntenic connection, we can
# add the exact scaffold name later, for example:
#
# flips <- data.frame(
#   genome = focal_genome,
#   chr = "contig_X"
# )
#
# and rerun plot_riparian() with:
#
# invertTheseChrs = flips
#
# This affects visualization/orientation only and does not alter the assembly.


# ============================================================================
# Final checks
# ============================================================================

expected_outputs <- c(
  pairwise_pdf,
  all_genomes_pdf
)

missing_outputs <- expected_outputs[
  !file.exists(expected_outputs)
]

if (length(missing_outputs) > 0) {

  stop(
    paste(
      "Missing expected plot(s):",
      paste(
        missing_outputs,
        collapse = ", "
      )
    )
  )

}


message("")
message("============================================")
message("Synteny plots completed successfully")
message("============================================")

message("")
message("Plot 09:")
message(pairwise_pdf)

message("")
message("Plot 10:")
message(all_genomes_pdf)

message("")
message("GENESPACE analysis copies:")
message(riparian_dir)

message("")
message(
  "Next step: inspect Plot 09 for reversed scaffolds or unusual syntenic connections."
)

