---
title: "Mitochondrial Genome Assembly"
author: "Tamoghno Paul"
date: '2026'
output:
  html_document:
    toc: true
    toc_float: true
    theme: cosmo
    highlight: tango
  pdf_document:
    toc: true
---

# Appendix

## A.1 SLURM Job Script — SRA Prefetch

**File:** `prefetch_job.sh`

```{bash eval=FALSE}
#!/bin/bash
#SBATCH --get-user-env
#SBATCH --partition=prevail
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00
#SBATCH -J fuckumean
prefetch -X 30G $1
```

---

### SBATCH Directives

| Directive | Value | Description |
|---|---|---|
| `--get-user-env` | — | Inherits the submitting user's environment |
| `--partition` | `prevail` | Target cluster partition |
| `--cpus-per-task` | `1` | Single CPU core per task |
| `--mem-per-cpu` | `4763 MB` | Memory per CPU (~4.65 GB) |
| `--time` | `48:00:00` | Maximum wall-clock time of 48 hours |
| `-J` | `fuckumean` | Job name |

---

### Usage

```{bash eval=FALSE}
sbatch prefetch_job.sh <SRR_ACCESSION>
```

**Example accessions used with this script:**

| Accession | Notes |
|---|---|
| `SRR38417560` | |
| `SRR30752828` | |
| `SRR36678954` | |
| `SRR24465296` | |

---

### Notes

- The `-X 30G` flag sets the maximum download size per accession to **30 GB**.
- Accessions are passed as a positional argument (`$1`) at submission time, allowing the same script to be reused across all samples.
- To submit all four accessions, run `sbatch` once per accession or wrap in a loop:

```{bash eval=FALSE}
for SRR in SRR38417560 SRR30752828 SRR36678954 SRR24465296; do
    sbatch prefetch_job.sh $SRR
done
```

## A.2 Mitochondrial Genome Assembly Pipeline

This pipeline maps paired-end reads to a reference mitochondrial sequence, extracts mapped reads, assembles them with SPAdes, and retrieves the mitogenome contig.

---

### A.2.1 Pipeline Overview

| Step | Tool | Purpose |
|---|---|---|
| 1 | `bwa index` | Index the reference mitochondrial FASTA |
| 2 | `bwa mem` + `samtools` | Map reads, filter unmapped, sort |
| 3 | `samtools sort -n` | Sort BAM by read name |
| 4 | `samtools fastq` | Extract mapped reads back to FASTQ |
| 5 | `grep -c` | Count recovered read pairs |
| 6 | `spades.py` | Assemble mitochondrial reads |
| 7 | `grep ">"` | Inspect contigs — **note the NODE_1 name here** |
| 8 | `samtools faidx` | Extract the mitogenome contig |
| 9 | `wc -c` | Verify genome length |

---

### A.2.2 Step-by-Step Commands

> ⚠️ **Before running each accession:** after SPAdes finishes, run the `grep ">"` command on `contigs.fasta` to identify the correct `NODE_1` name — it **will differ** between accessions.

---

#### SRR24465296

```{bash eval=FALSE}
# Step 1 — Index reference (only needs to be done once)
bwa index sequence2.fasta

# Step 2 — Map, filter unmapped reads, sort
bwa mem -t 8 sequence2.fasta \
    ~/mito_project/SRR24465296_1.fastq \
    ~/mito_project/SRR24465296_2.fastq \
    | samtools view -b -F 4 \
    | samtools sort > SRR24465296_mt.bam

# Step 3 — Sort by read name
samtools sort -n SRR24465296_mt.bam -o SRR24465296_mt_namesorted.bam

# Step 4 — Extract mapped reads to FASTQ
samtools fastq \
    -1 SRR24465296_mt_R1.fastq \
    -2 SRR24465296_mt_R2.fastq \
    -0 /dev/null -s /dev/null \
    SRR24465296_mt_namesorted.bam

# Step 5 — Count read pairs
grep -c "^@" SRR24465296_mt_R1.fastq

# Step 6 — Assemble
spades.py --careful \
    -1 SRR24465296_mt_R1.fastq \
    -2 SRR24465296_mt_R2.fastq \
    -o SRR24465296_spades_out

# Step 7 — CHECK CONTIG NAMES HERE before proceeding
grep ">" SRR24465296_spades_out/contigs.fasta

# Step 8 — Extract mitogenome contig (update NODE name from Step 7)
samtools faidx SRR24465296_spades_out/contigs.fasta \
    NODE_1_length_16892_cov_413.120369 > SRR24465296_mitogenome.fasta

# Step 9 — Verify
grep ">" SRR24465296_mitogenome.fasta
cat SRR24465296_mitogenome.fasta | grep -v ">" | tr -d '\n' | wc -c
```

---

#### SRR38417560

```{bash eval=FALSE}
bwa mem -t 8 sequence2.fasta \
    ~/mito_project/SRR38417560_1.fastq \
    ~/mito_project/SRR38417560_2.fastq \
    | samtools view -b -F 4 \
    | samtools sort > SRR38417560_mt.bam

samtools sort -n SRR38417560_mt.bam -o SRR38417560_mt_namesorted.bam

samtools fastq \
    -1 SRR38417560_mt_R1.fastq \
    -2 SRR38417560_mt_R2.fastq \
    -0 /dev/null -s /dev/null \
    SRR38417560_mt_namesorted.bam

grep -c "^@" SRR38417560_mt_R1.fastq

spades.py --careful \
    -1 SRR38417560_mt_R1.fastq \
    -2 SRR38417560_mt_R2.fastq \
    -o SRR38417560_spades_out

# CHECK CONTIG NAMES HERE before proceeding
grep ">" SRR38417560_spades_out/contigs.fasta

# Update NODE name below based on grep output above
samtools faidx SRR38417560_spades_out/contigs.fasta \
    NODE_?_length_?_cov_? > SRR38417560_mitogenome.fasta

grep ">" SRR38417560_mitogenome.fasta
cat SRR38417560_mitogenome.fasta | grep -v ">" | tr -d '\n' | wc -c
```

---

#### SRR30752828

```{bash eval=FALSE}
bwa mem -t 8 sequence2.fasta \
    ~/mito_project/SRR30752828_1.fastq \
    ~/mito_project/SRR30752828_2.fastq \
    | samtools view -b -F 4 \
    | samtools sort > SRR30752828_mt.bam

samtools sort -n SRR30752828_mt.bam -o SRR30752828_mt_namesorted.bam

samtools fastq \
    -1 SRR30752828_mt_R1.fastq \
    -2 SRR30752828_mt_R2.fastq \
    -0 /dev/null -s /dev/null \
    SRR30752828_mt_namesorted.bam

grep -c "^@" SRR30752828_mt_R1.fastq

spades.py --careful \
    -1 SRR30752828_mt_R1.fastq \
    -2 SRR30752828_mt_R2.fastq \
    -o SRR30752828_spades_out

# CHECK CONTIG NAMES HERE before proceeding
grep ">" SRR30752828_spades_out/contigs.fasta

# Update NODE name below based on grep output above
samtools faidx SRR30752828_spades_out/contigs.fasta \
    NODE_?_length_?_cov_? > SRR30752828_mitogenome.fasta

grep ">" SRR30752828_mitogenome.fasta
cat SRR30752828_mitogenome.fasta | grep -v ">" | tr -d '\n' | wc -c
```

---

#### SRR36678954

```{bash eval=FALSE}
bwa mem -t 8 sequence2.fasta \
    ~/mito_project/SRR36678954_1.fastq \
    ~/mito_project/SRR36678954_2.fastq \
    | samtools view -b -F 4 \
    | samtools sort > SRR36678954_mt.bam

samtools sort -n SRR36678954_mt.bam -o SRR36678954_mt_namesorted.bam

samtools fastq \
    -1 SRR36678954_mt_R1.fastq \
    -2 SRR36678954_mt_R2.fastq \
    -0 /dev/null -s /dev/null \
    SRR36678954_mt_namesorted.bam

grep -c "^@" SRR36678954_mt_R1.fastq

spades.py --careful \
    -1 SRR36678954_mt_R1.fastq \
    -2 SRR36678954_mt_R2.fastq \
    -o SRR36678954_spades_out

# CHECK CONTIG NAMES HERE before proceeding
grep ">" SRR36678954_spades_out/contigs.fasta

# Update NODE name below based on grep output above
samtools faidx SRR36678954_spades_out/contigs.fasta \
    NODE_?_length_?_cov_? > SRR36678954_mitogenome.fasta

grep ">" SRR36678954_mitogenome.fasta
cat SRR36678954_mitogenome.fasta | grep -v ">" | tr -d '\n' | wc -c
```

---

### A.2.3 NODE Selection Guide

After running `grep ">" <accession>_spades_out/contigs.fasta`, you will see output like:

```
>NODE_1_length_16892_cov_413.120369
>NODE_2_length_312_cov_5.210000
>NODE_3_length_98_cov_1.430000
```

Select the node that matches expected mitogenome criteria:

| Criterion | Expected value |
|---|---|
| **Length** | ~15,000 – 17,000 bp |
| **Coverage** | Significantly higher than other nodes |
| **Node number** | Usually `NODE_1` |

## A.3 Alignment, Cleaning, and Phylogenetic Inference

---

### A.3.1 Pipeline Overview

| Step | Tool | Purpose |
|---|---|---|
| 1 | `cat` | Concatenate all mitogenomes into one FASTA |
| 2 | `mafft` | Multiple sequence alignment |
| 3 | `pxclsq` | Clean alignment at multiple stringency thresholds |
| 4 | `iqtree3` | Phylogenetic inference for each cleaned alignment |

---

### A.3.2 Step-by-Step Commands

#### Step 1 — Concatenate Mitogenomes

```{bash eval=FALSE}
cat SRR30752828_mitogenome.fasta \
    SRR24465296_mitogenome.fasta \
    SRR18911050_mitogenome.fasta \
    SRR38417560_mitogenome.fasta \
    > all_mitogenome.fasta
```

---

#### Step 2 — Multiple Sequence Alignment with MAFFT

```{bash eval=FALSE}
mafft --auto --adjustdirection all_mitogenome.fasta > all_mitogenome_aligned.aln
```

| Flag | Description |
|---|---|
| `--auto` | Automatically selects the best alignment strategy based on data size |
| `--adjustdirection` | Automatically reverse-complements sequences if needed |

---

#### Step 3 — Alignment Cleaning with pxclsq

Three thresholds are tested to assess the effect of missing data filtering on the phylogeny:

```{bash eval=FALSE}
pxclsq -s all_mitogenomes_aligned.aln -o all_mitogenomes_cleaned_03.aln-cln -p 0.3
pxclsq -s all_mitogenomes_aligned.aln -o all_mitogenomes_cleaned_04.aln-cln -p 0.4
pxclsq -s all_mitogenomes_aligned.aln -o all_mitogenomes_cleaned_05.aln-cln -p 0.5
```

| Output file | `-p` threshold | Meaning |
|---|---|---|
| `all_mitogenomes_cleaned_03.aln-cln` | `0.3` | Sites with ≥30% occupancy retained |
| `all_mitogenomes_cleaned_04.aln-cln` | `0.4` | Sites with ≥40% occupancy retained |
| `all_mitogenomes_cleaned_05.aln-cln` | `0.5` | Sites with ≥50% occupancy retained |

> ⚠️ Higher `-p` values = stricter cleaning = fewer sites retained but less missing data.

---

#### Step 4 — Phylogenetic Inference with IQ-TREE 3

One tree is inferred per cleaned alignment:

```{bash eval=FALSE}
iqtree3 -m TEST --merge \
    -s all_mitogenomes_cleaned_03.aln-cln \
    -T AUTO -B 1000 -wbtl \
    --seqtype DNA \
    --prefix IQtree3_chl_alignments_phyx03

iqtree3 -m TEST --merge \
    -s all_mitogenomes_cleaned_04.aln-cln \
    -T AUTO -B 1000 -wbtl \
    --seqtype DNA \
    --prefix IQtree3_chl_alignments_phyx04

iqtree3 -m TEST --merge \
    -s all_mitogenomes_cleaned_05.aln-cln \
    -T AUTO -B 1000 -wbtl \
    --seqtype DNA \
    --prefix IQtree3_chl_alignments_phyx05
```

| Flag | Description |
|---|---|
| `-m TEST` | Automatically selects the best substitution model |
| `--merge` | Tests for merging partitions with similar rates |
| `-T AUTO` | Automatically determines optimal number of threads |
| `-B 1000` | Ultrafast bootstrap with 1000 replicates |
| `-wbtl` | Writes bootstrap trees with branch lengths |
| `--seqtype DNA` | Specifies DNA input |
| `--prefix` | Output file prefix for each run |

---

### A.3.3 Output Files per IQ-TREE Run

Each `--prefix` generates the following key files:

| Extension | Contents |
|---|---|
| `.treefile` | Best-scoring ML tree |
| `.iqtree` | Full run report including selected model |
| `.log` | Run log |
| `.contree` | Consensus tree with bootstrap values |
| `.ufboot` | Bootstrap trees with branch lengths |
