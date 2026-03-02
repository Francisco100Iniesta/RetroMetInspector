# RetroMetInspector

RetroMetInspector is a multi-omics extension of the RetroInspector framework (https://github.com/javiercguard/retroinspector) for the functional characterization of transposable element (TE) insertions using long-read sequencing data.

This workflow integrates TE detection with haplotype-aware methylation profiling, enabling allele-specific epigenetic interpretation of insertion events.

---

## 🔬 Conceptual Overview

Transposable elements (TEs) are major contributors to genome variability. While structural variant-based tools can detect TE insertions, their potential epigenetic impact remains largely unexplored.

RetroMetInspector bridges structural genomics and epigenomics by combining:

- TE insertion detection  
- Haplotype assignment  
- Allele-specific methylation quantification  
- Differential methylation analysis at insertion loci  

This enables functional interpretation of heterozygous insertions.

---

## 🧬 Functional Extension Overview

RetroMetInspector is a downstream analytical framework built on top of RetroInspector TE calls.

Rather than performing de novo TE detection, RetroMetInspector:

1. Processes TE insertions previously identified by RetroInspector
2. Performs haplotype-aware characterization of insertion loci
3. Computes allele-specific methylation profiles
4. Generates multi-layer methylation visualizations
5. Enables inter-sample comparison of insertion-associated methylation states

---

## 🧱 Relationship to RetroInspector

TE insertion discovery is performed using the RetroInspector workflow:
https://github.com/javiercguard/retroinspector

RetroMetInspector takes these reported insertion sites as input and performs extended epigenomic characterization.

---

## 🔍 Core Functionalities

### 1️⃣ Haplotype-aware insertion genotyping

- Assigns insertions to specific haplotypes
- Determines the insertion-bearing allele
- Refines genotyping accuracy at heterozygous loci
- Improves allele-specific interpretation
<img width="644" height="506" alt="image" src="https://github.com/user-attachments/assets/00b58ac8-e326-4d65-8499-17c5664b6b16" />

---

### 2️⃣ Customizable methylation window analysis

For each insertion locus, the tool:

- Extracts user-defined genomic windows surrounding the insertion
- Computes methylation levels per haplotype
- Detects allele-specific methylation differences
- Applies configurable filtering thresholds

This allows detection of subtle methylation shifts (e.g., MINmetDiff).
<img width="1768" height="178" alt="image" src="https://github.com/user-attachments/assets/561f34cf-b933-4d12-9ca3-975aea589e15" />
<img width="1760" height="121" alt="image" src="https://github.com/user-attachments/assets/8ffe50cf-05e9-4ce6-8f44-ca0e11af440d" />


---

### 3️⃣ Visualization via Methylartist

Methylation plots are generated using:
https://github.com/adamewing/methylartist

The workflow automatically:

- Produces methylation profiles across selected windows
- Generates haplotype-separated plots
- Outputs publication-ready figures
- Integrates plots into interactive HTML reports
<img width="870" height="449" alt="image" src="https://github.com/user-attachments/assets/2593d1bb-1a6b-493a-8161-1d81450c93c2" />
---

### 4️⃣ Custom contig reconstruction of insertion loci

To directly assess methylation within the inserted element:

- Custom contigs are constructed by incorporating the insertion sequence into the reference locus
- Reads are re-evaluated against these reconstructed regions
- Methylation within the inserted TE itself is computed and plotted

This enables direct epigenetic profiling of the insertion sequence.
<img width="1051" height="566" alt="image" src="https://github.com/user-attachments/assets/4be72c02-1c92-4842-9f8d-f86422bb46fd" />

---

### 5️⃣ Cross-sample comparison mode ("compare" mode)

RetroMetInspector includes a comparative module that:

- Identifies shared or unique insertions between samples
- Compares methylation states between:
    - Samples carrying the insertion
    - Samples lacking the insertion
- Reports allele-specific methylation differences across conditions

This enables functional interpretation of insertion-associated epigenetic variation across individuals or experimental conditions.
![Sin título-1](https://github.com/user-attachments/assets/0a38dafd-3caf-430e-bec9-572701e09050)

---

## 🧠 Biological Interpretation Enabled

RetroMetInspector enables questions such as:

- Does a heterozygous TE insertion alter local methylation only in the insertion-bearing allele?
- Is the inserted element itself methylated?
- Are methylation differences preserved across individuals?
- Does insertion presence correlate with allele-specific epigenetic shifts?

---
## Input requirements

RetroMetInspector runs **downstream of RetroInspector** and requires haplotype-resolved long-read data with methylation information.

**Mandatory inputs**
- **Reference genome** (`.fa/.fasta`) and **annotation** (`.gtf/.gff`)
- **Aligned long-read BAM(s)** for each sample

**BAM must contain**
- **Methylation tags**: `MM` and `ML` (typically produced by ONT **Dorado** basecaller)
- **Phasing / haplotype tag**: `HP` with reads assigned to **HP1 / HP2** (e.g., generated with **WhatsHap**)

Without `MM/ML` the workflow cannot quantify methylation, and without `HP` it cannot perform allele-specific (haplotype-aware) analysis.

**Sample naming (config)**
- The `samples` list must match your BAM identifiers (no spaces recommended).
- In **compare mode**, each entry in `comparisons` must contain **exactly two** sample names.

```yaml
samples: { sampleA, sampleB, sampleC }

comparisons:
  [
    ['sampleA', 'sampleB']
  ]


## Launch

RetroMetInspector is executed using **Snakemake** with Conda/Mamba environments.

It is recommended to run the workflow on an HPC system.

### Environment activation

```{bash}
source ~/.bashrc
mamba activate snakemake
snakemake --use-conda \
  --conda-frontend mamba \
  -p \
  --configfile config.default.yaml \
  -c 46




referenceGenome: "data/hg38.fa"
humanGTF: "data/gencode.v49.annotation.sorted.gtf.gz"
species: "homo_sapiens"

bam_directory: "alns"
outputPath: "."
allPrefix: "project_prefix"

# --- Methylation analysis parameters ---
MetWindows: 1000
MINmetDiff: 0.25
MinreadsSupport: 5

threads: 32
mode: "full"

samples:
  { sampleA, sampleB }

comparisons:
  [
    ['sampleA', 'sampleB']
  ]
