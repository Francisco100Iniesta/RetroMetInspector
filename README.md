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

---

### 2️⃣ Customizable methylation window analysis

For each insertion locus, the tool:

- Extracts user-defined genomic windows surrounding the insertion
- Computes methylation levels per haplotype
- Detects allele-specific methylation differences
- Applies configurable filtering thresholds

This allows detection of subtle methylation shifts (e.g., MINmetDiff).

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
<img width="870" height="449" alt="image" src="https://github.com/user-attachments/assets/2593d1bb-1a6b-493a-8161-1d81450c93c2" />
<img width="870" height="449" alt="image" src="https://github.com/user-attachments/assets/2593d1bb-1a6b-493a-8161-1d81450c93c2" />

---

### 4️⃣ Custom contig reconstruction of insertion loci

To directly assess methylation within the inserted element:

- Custom contigs are constructed by incorporating the insertion sequence into the reference locus
- Reads are re-evaluated against these reconstructed regions
- Methylation within the inserted TE itself is computed and plotted

This enables direct epigenetic profiling of the insertion sequence.

---

### 5️⃣ Cross-sample comparison mode ("compare" mode)

RetroMetInspector includes a comparative module that:

- Identifies shared or unique insertions between samples
- Compares methylation states between:
    - Samples carrying the insertion
    - Samples lacking the insertion
- Reports allele-specific methylation differences across conditions

This enables functional interpretation of insertion-associated epigenetic variation across individuals or experimental conditions.

---

## 🧠 Biological Interpretation Enabled

RetroMetInspector enables questions such as:

- Does a heterozygous TE insertion alter local methylation only in the insertion-bearing allele?
- Is the inserted element itself methylated?
- Are methylation differences preserved across individuals?
- Does insertion presence correlate with allele-specific epigenetic shifts?

---
## ⚙️ Configuration Example

```yaml
referenceGenome: "data/hg38.fa"
humanGTF: "data/yourGTF"
species: "homo_sapiens"
bam_directory: "alns"
outputPath: "."
allPrefix: "project_prefix"

# Methylation analysis parameters
MetWindows: 1000
MINmetDiff: 0.25
MinreadsSupport: 5
