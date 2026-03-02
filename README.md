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

## 🧱 Relationship to RetroInspector

RetroMetInspector builds upon the TE detection framework implemented in RetroInspector:

https://github.com/javiercguard/retroinspector

RetroInspector performs:

- Detection of TE insertions and deletions
- Structural variant annotation
- Configurable filtering by read support
- Comparative analysis between samples

RetroMetInspector extends this framework by adding:

- Haplotype-resolved insertion characterization  
- Allele-specific methylation quantification  
- Detection of methylation differences between haplotypes  
- Customizable methylation windows around insertion sites  
- Advanced visualization of methylation profiles  
- Functional multi-omic interpretation  

TE detection itself is performed using the RetroInspector-based workflow, while RetroMetInspector introduces additional epigenetic and haplotype-aware modules.

---

## 🧬 Key Features

✔ TE insertion detection (via RetroInspector framework)  
✔ Haplotype assignment of insertions  
✔ Allele-specific methylation quantification  
✔ Differential methylation analysis (Haplotype 1 vs Haplotype 2)  
✔ Customizable genomic windows for methylation analysis  
✔ Configurable filtering thresholds  
✔ Automated HTML report generation  
✔ Publication-ready methylation plots  

This enables biological questions such as:

> Can a heterozygous TE insertion alter the methylation state of the surrounding genomic region specifically in the insertion-bearing allele?

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
