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
- **IGV coordinate file generation**: After pipeline completion, the workflow generates a coordinate file containing all detected insertion loci. This file enables direct visualization in IGV and facilitates systematic screenshot capture of insertion events. Detailed instructions for IGV-based inspection and screenshot generation are provided in a later section.
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

### Mandatory inputs

- **Reference genome** (`.fa/.fasta`)
- **Gene annotation** (`.gtf/.gff`)
- **Aligned long-read BAM files** for each sample
- **TE insertion coordinates identified by RetroInspector**

### BAM requirements

Input BAM files must contain:

- **Methylation tags**: `MM` and `ML`  
  (typically generated by ONT *Dorado* basecaller)

- **Haplotype tag**: `HP`  
  - `HP:1` → Haplotype 1  
  - `HP:2` → Haplotype 2  
  (typically assigned using *WhatsHap*)

Without `MM/ML` tags, methylation quantification cannot be performed.  
Without `HP` tags, allele-specific methylation analysis cannot be computed.

### Sample naming (configuration file)

The `samples` list must match your BAM identifiers (no spaces recommended).

In **compare mode**, each entry in `comparisons` must contain exactly two samples.

```{yaml}
samples: { sampleA, sampleB, sampleC, ... }

comparisons:
  [
    ['sampleA', 'sampleB']
  ]
```

### Clone repository
```{yaml}
git clone https://github.com/Francisco100Iniesta/RetroMetInspector
```
### Launch pipeline

```{yaml}
source ~/.bashrc
mamba activate snakemake
snakemake --use-conda \
  --conda-frontend mamba \
  -p \
  --configfile config.default.yaml \
  -c 46
```

### Config File example
```{yaml}
referenceGenome: "data/hg38.fa"
humanGTF: "data/gencode.v49.annotation.sorted.gtf.gz"
species: "homo_sapiens"
bam_directory: "alns"
outputPath: "." # Dont't include the trailing "/"
allPrefix: "organoid"
#params for retromet analysis
MetWindows: 1000
MINmetDiff: 0.25
MinreadsSupport: 5
# Parameters
threads: 32
mode: "full"
minimumReadSupport: 1
insertionDistanceLimitIntraPatient: 200
insertionDistanceLimitInterPatient: 200
survivorInsertionDistanceLimitIntraPatient: 200
survivorInsertionDistanceLimitInterPatient: 500
enrichmentSignificanceThreshold: 0.05
callers: ["cuteSV", "sniffles2"]
datasets: []

# Keep certain files
keepRds: True

samples: {SampleA, SampleB, SampleC}

comparisons:
  [
   ['SampleA', 'SampleB'],
  ]

```
### Parameter description

- **referenceGenome**: Path to the reference genome FASTA file used as genomic coordinate system.

- **humanGTF**: Gene annotation file (GTF/GFF) used for feature annotation.

- **species**: Species identifier (e.g., *homo_sapiens*).

- **bam_directory**: Directory containing phased long-read BAM files with MM/ML (methylation) and HP (haplotype) tags.

- **outputPath**: Directory where results will be written (do not include a trailing slash).

- **allPrefix**: Prefix used for naming output files.

- **MetWindows**: Size (bp) of the methylation window analyzed on each side of the insertion breakpoint.

- **MINmetDiff**: Minimum methylation difference (0–1 scale) between haplotypes required to report differential methylation.

- **MinreadsSupport**: Minimum number of supporting reads per haplotype required to compute and plot methylation.

- **threads**: Number of computational threads available.

- **mode**: Workflow execution mode (e.g., `"full"` for complete haplotype-aware methylation analysis).

- **minimumReadSupport**: Minimum read support required for TE insertion calls.

- **insertionDistanceLimitIntraPatient**: Distance threshold (bp) to merge insertions within the same sample.

- **insertionDistanceLimitInterPatient**: Distance threshold (bp) to merge insertions between samples.

- **survivorInsertionDistanceLimitIntraPatient**: SURVIVOR merging threshold within sample.

- **survivorInsertionDistanceLimitInterPatient**: SURVIVOR merging threshold between samples.

- **enrichmentSignificanceThreshold**: Statistical significance cutoff for enrichment analysis.

- **callers**: Structural variant callers used for TE detection (e.g., cuteSV, sniffles2).

- **datasets**: Optional grouping variable for multi-dataset analysis.

- **keepRds**: If `True`, intermediate R objects are preserved for downstream reuse.

- **samples**: List of sample identifiers (must match BAM file names).

- **comparisons**: Pairwise sample comparisons used in compare mode (exactly two samples per comparison).
  
## IGV Screenshot generation
A script is available in the IGV DIRECTORY to generate screenshots of all haplotyped transposition events in a comparative way.

The script can be runned after the pipeline run.
```{yaml}
./screenshot_igv.sh {bat file generated in the IGV directory} {File path to BAM 1} {File path to BAM 2} {Output path to screenshot directory (blank by default)}
```
Screenshot output: Generates a directory containing one image per insertion event, where each screenshot file is named using the corresponding insertion ID.

<img width="443" height="345" alt="image" src="https://github.com/user-attachments/assets/2c7c12e0-ff63-48be-b609-e988c8c6fde0" />

