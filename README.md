# Single-cell WGS allele-specific copy-number analysis with ASCAT.sc

Scripts for PicoPLEX single-cell whole-genome sequencing preprocessing, allele counting, Beagle phasing, ASCAT.sc copy-number analysis and downstream visualisation.

## Workflow

| Stage | Script | Description |
| --- | --- | --- |
| 1 | `01_preprocessing/` | FASTQ trimming, GRCh38 alignment and duplicate marking |
| 2 | `02_allele_counts/` | Bulk and single-cell allele counting |
| 3 | `02_phasing/` | Confident-SNP bulk phasing with Beagle |
| 4 | `03_ascat/` | Allele-specific copy-number analysis with ASCAT.sc |
| 5 | `04_post_ascat/` | Copy-number and BAF visualisation |

## Usage

Configure the reference files, input directories and output paths in `config/example_environment.sh`.

```bash
source config/example_environment.sh
```

### Preprocessing

```bash
bash 01_preprocessing/01_trim_PicoPLEX.sh /path/to/raw_fastqs
qsub -V 01_preprocessing/02_align_GRCh38.sh
qsub -V 01_preprocessing/03_mark_duplicates.sh
```

### Allele counting

Bulk and single-cell allele counting are submitted independently using:

```bash
qsub -V 02_allele_counts/submit_bulk_allele_counts_SGE.sh
qsub -V 02_allele_counts/submit_single_cell_allele_counts_SGE.sh
```

### Beagle phasing

```bash
qsub -V -t 1-23 02_phasing/submit_beagle_SGE.sh
```

### ASCAT.sc

```bash
qsub -V -t 1-4 03_ascat/submit_ASCATsc_SGE.sh
```

The pipeline uses bulk Beagle-phased SNPs and independently generated single-cell allele counts for ASCAT.sc analysis.

## Requirements

- R and Python 3
- ASCAT.sc and associated R packages
- BWA-MEM, SAMtools, Picard, fastp and seqkit
- alleleCounter and Beagle
- Apptainer and Java
- GRCh38 reference genome and SNP resources
- SGE/Grid Engine for cluster submission

Reference datasets and sequencing files are not included in this repository.