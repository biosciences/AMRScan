# AMRScan 🧬

**AMRScan** is a hybrid bioinformatics toolkit written in **R** and **Nextflow** for the rapid detection of antimicrobial resistance (AMR) genes from next-generation sequencing (NGS) data. It aligns sequencing reads against a curated AMR gene database (e.g., CARD) using BLAST and reports likely resistance determinants. It also supports both single-sample scripting for rapid testing and scalable, containerized pipelines for high-throughput environments.

---

## 📄 Project Links
- 💻 [View Source Code on GitHub](https://github.com/biosciences/AMRScan): Explore the full repository
- 🧬 [Live Report (GitHub Pages)](https://biosciences.github.io/AMRScan/index.html): View the interactive HTML output

---

## 🚀 Features

- Antimicrobial resistance gene detection using BLAST
- Compatible with reference databases such as CARD
- Two modes of use:
  - **R script**: Simple, linear script for single samples
  - **Nextflow workflow**: Scalable, reproducible, and parallelizable pipeline
- Supports FASTQ input, automatic conversion, quality control, and output summarization
- Read quality control from FASTQ inputs
- Conversion to FASTA and BLAST-compatible formats
- Automated `blastx` comparison to AMR gene databases
- Top-hit summarization of resistance genes per read
- Output in tidy `.csv` format for downstream analysis

---

## 🤔 When to Use Which Version?

| Scenario                                             | Use R Script (`AMRScan_standalone.R`) | Use Nextflow Workflow (`main.nf`)  |
|------------------------------------------------------|----------------------------------------|-------------------------------------|
| Small dataset / one or two samples                   | ✅ Yes                                 | ✅ Yes                              |
| Educational / demonstration setting                  | ✅ Yes                                 | ❌ No                               |
| Quick local prototyping with minimal dependencies    | ✅ Yes                                 | ❌ No                               |
| Batch processing of many samples                     | ❌ No                                  | ✅ Yes                              |
| HPC, cloud, or Docker/Conda-based reproducibility    | ❌ No                                  | ✅ Yes                              |
| Automated multi-step workflows (BLAST, parsing, etc) | ❌ No                                  | ✅ Yes                              |

## ✅ How the files are used in the dual-mode setup

| File              | Role                                                    | Executed in R? | Executed by Nextflow? |
|-------------------|---------------------------------------------------------|----------------|-----------------------|
| `AMRScan.R`       | Full R script version (standalone)                      | ✅ Yes         | ❌ No                 |
| `AMRScan.nf`      | Nextflow pipeline script (main workflow logic)          | ❌ No          | ✅ Yes                |
| `nextflow.config` | Configuration file for resource settings and parameters | ❌ No          | ✅ Yes                |

---
## 🚀 Usage

You can run AMRScan in two modes:

### 🧪 Option 1: R Script

#### Requirements
- [R (>= 4.0)](https://cran.r-project.org/)
- R packages: `dplyr`, `RCurl`
- Bioconductor packages: `ShortRead`, `Biostrings`
- [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

#### 📥 Installation

Ensure the above requried R, R packages and NCBI BLAST+ are installed.

Clone the repo:
```bash
git clone https://github.com/biosciences/AMRScanR.git
cd AMRScanR
```

Install dependencies in R:
```r
install.packages(c("dplyr", "RCurl"))
BiocManager::install(c("ShortRead", "Biostrings"))
```

Run
```r
source("scripts/AMRScan_standalone.R")
```
The script will:
•Clean your FASTQ input
•Convert reads to FASTA
•BLAST against the CARD database (auto-downloads and formats if needed)
•Generate a AMR_hits_summary.csv in the results/ folder

#### Output
- Results written to `results/AMR_hits_summary.csv`

### 🧬 Option 2: Nextflow Workflow

#### Requirements
- Nextflow (>= 22.04.0)
- Docker or Conda (recommended for reproducibility)
- Tools: `wget`, `blastx`, `seqtk`, `Rscript`

#### Run

```bash
# Basic usage
nextflow run workflow/AMRScan.nf --input data/GCF_037966445.1_ASM3796644v1_genomic.fna --outdir results

# With custom parameters
nextflow run workflow/AMRScan.nf \
    --input data/*.fna \
    --outdir results \
    --evalue 1e-10 \
    --threads 8 \
    --min_quality 25
```

#### Configuration
`workflow/nextflow.config` includes basic container settings. You can customize profiles for HPC/cloud usage.

---

## 📂 Folder Structure

```
AMRScan/
├── AMRScan_Nextflow.Rmd        # Nextflow Markdown for live demo
├── AMRScan_R.Rmd               # R Markdown for live demo
├── data/                       # Input FASTA file(s)
│   └── GCF_037966445.1_ASM3796644v1_genomic.fna
├── db/                    # BLAST database folder
│   └── protein_fasta_protein_homolog_model.fasta
├── DESCRIPTION                 # Package or metadata description
├── docs/                       # GitHub Pages rendered site
│   ├── AMRScan_Nextflow.html
│   ├── AMRScan_R.html
│   └── index.html
├── LICENSE                     # Software license
├── paper.md                    # JOSS submission manuscript
├── paper.bib                   # Bibliography for JOSS paper
├── README.md                   # Project documentation
├── results/                    # Output results directory
├── blast/                      # Intermediate BLAST results
│   └── GCF_037966445.1_ASM3796644v1_genomic_blast_results.tsv
├── final/                      # Final processed results
│   ├── GCF_037966445.1_ASM3796644v1_genomic_AMR_hits_detailed.csv
│   └── GCF_037966445.1_ASM3796644v1_genomic_AMR_hits_summary.csv
├── preprocessed/               # Preprocessed intermediate files
│   ├── GCF_037966445.1_ASM3796644v1_genomic.fasta
│   └── GCF_037966445.1_ASM3796644v1_genomic_preprocess.log
├── scripts/                    # Scripts for analysis and setup
│   └── AMRScan.R               # Standalone R version
├── tests/                      # Test files
│   ├── test-run.R              # R test resource
├── workflow/                   # Nextflow pipeline and configuration
│   ├── AMRScan.nf              # Nextflow script
│   └── nextflow.config         # Configuration file
```

---

## 📂 Input Example

Download fasta data from:

```
https://www.ncbi.nlm.nih.gov/assembly/GCA_037966445.1
```

The fasta data source from:

```
Munim, M. A., Tanni, A. A., Hossain, M. M., Chakma, K., Mannan, A., Islam, S. M. R., Tiwari, J. G., & Gupta, S. D. (2024). *Whole genome sequencing of multidrug-resistant Klebsiella pneumoniae from poultry in Noakhali, Bangladesh: Assessing risk of transmission to humans in a pilot study*. Comparative Immunology, Microbiology and Infectious Diseases, 114, 102246. https://doi.org/10.1016/j.cimid.2024.102246
```
Download the comprehensive antibiotic resistance database from:
https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2

Uncompress above tar.bz2 file, copy protein_fasta_protein_homolog_model.fasta to data folder


## 📊 Example Output

| Query                | Subject     | Identity | Length | Evalue   | Bitscore | Annotation                    |
|----------------------|-------------|----------|--------|----------|----------|-------------------------------|
| NZ_JBBPBW010000028.1 | ARO:3003039 | 40.839   | 453    | 1.55e-72 | 252      | OprA [Pseudomonas aeruginosa] |
| NZ_JBBPBW010000035.1 | ARO:3004826 | 100.000  | 285    | 0        | 587      | LAP-2 [Enterobacter cloacae]  |

---

## 📊 Live Report

An interactive demonstration of AMRScan can be viewed here:

👉 [View the AMRScan Demo Report](https://your-username.github.io/AMRScan/)

This report is generated from `AMRScan.Rmd` and showcases a complete run-through using example input files, result parsing, and summary.

---

## 📊 Test

Test AMRScan.R
```r
source("tests/test-run.R")
```

---

## 🧾 Citation

If you use **AMRScan** in your research, please cite the associated JOSS paper (under review):

```
Lai, K. (2025). AMRScan: A hybrid R and Nextflow toolkit for rapid antimicrobial resistance gene detection. Journal of Open Source Software. (under review)
```

---

## 🪪 License

MIT © 2025 Kaitao Lai
