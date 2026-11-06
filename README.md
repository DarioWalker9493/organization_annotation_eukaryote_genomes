# Organization and Annotation of Eukaryote Genomes

Scripts, selected results and notes from the MSc course.

## Repository structure

- `scripts/` — analysis scripts
- `plots/` — selected figures
- `data/` — local input data and links to course files
- `results/` — analysis outputs
- `logs/` — Slurm output and error logs

## Week 1 — Transposable Element Annotation and Classification

### Genome assembly

The Flye assembly of the Lu-1 accession was used for transposable element annotation.

- Assembly size: **138,067,759 bp**
- Number of sequences: **99**
- TE annotation performed with **EDTA**
- Full-length LTR retrotransposons further classified with **TEsorter / REXdb-plant**

### EDTA annotation summary

EDTA annotated a total of **19,850,279 bp** as repetitive sequence, corresponding to **14.38% of the genome**.

| TE class / superfamily | Count | bp masked | Genome masked |
|---|---:|---:|---:|
| Gypsy | 5,232 | 7,076,529 | 5.13% |
| Helitron | 8,813 | 5,100,757 | 3.69% |
| Mutator | 2,488 | 1,799,501 | 1.30% |
| Copia | 1,557 | 1,754,341 | 1.27% |
| CACTA | 3,166 | 1,612,430 | 1.17% |
| Repeat fragments | 2,774 | 815,010 | 0.59% |
| L1 | 1,012 | 610,244 | 0.44% |
| LTR unknown | 715 | 449,893 | 0.33% |
| PIF/Harbinger | 722 | 297,560 | 0.22% |
| hAT | 419 | 195,365 | 0.14% |
| Tc1/Mariner | 55 | 69,590 | 0.05% |
| Jockey | 85 | 30,543 | 0.02% |
| Endogenous retrovirus | 40 | 28,981 | 0.02% |
| L2 | 18 | 8,591 | 0.01% |

**Main observation:** Gypsy elements occupy the largest fraction of the genome (**5.13%**), while Helitrons have the highest number of annotated copies (**8,813**).

The total number of interspersed repeat annotations was **27,099**, covering **14.38%** of the assembly.

### Full-length LTR classification

Full-length LTR retrotransposons identified by EDTA were classified with TEsorter using the `rexdb-plant` database.

The LTR identity distributions are strongly shifted toward high sequence identity, particularly around **0.95–1.00**, suggesting that many intact LTR elements represent relatively recent insertions. A smaller number of elements show identities around **0.80–0.90**, consistent with comparatively older insertions.

Different Copia and Gypsy clades show distinct identity distributions. Examples include strong high-identity signals in the Copia **Ale** and **Ivana** clades and the Gypsy **Athila** clade.

### Full-length LTR identity distributions

The identity distributions of full-length LTR retrotransposons are strongly shifted toward high identity values, particularly between 0.95 and 1.00. This suggests that many intact LTR elements represent relatively recent insertions.

#### Copia and Gypsy clade-level distributions

![Full-length LTR identity by Copia and Gypsy clade](plots/01_LTR_Copia_Gypsy_cladelevel.png)

Different clades show distinct identity distributions. Strong high-identity signals are visible, for example, for the Copia Ale and Ivana clades and the Gypsy Athila clade.

#### Global Copia versus Gypsy distribution

![Global Copia and Gypsy LTR identity distribution](plots/02_LTR_Copia_Gypsy_global.png)

Both Copia and Gypsy show a strong enrichment toward high LTR identity, while a smaller number of elements occur at lower identity values, consistent with comparatively older insertions.

### Workflow scripts

- `scripts/Week01_Script01_run_EDTA.sh` — EDTA genome-wide TE annotation
- `scripts/Week01_Script02_run_TEsorter_full_length_LTRs.sh` — classification of intact LTR retrotransposons
- `scripts/Week01_Script03_full_length_LTR_identity.R` — analysis and visualization of LTR identity
- `scripts/Week01_Script03_run_full_length_LTR_identity.sh` — Slurm wrapper for the R analysis