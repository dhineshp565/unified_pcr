# ONT Target Sequencing Pipeline - Usage Guide

## Table of Contents
1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Pipeline Steps](#pipeline-steps)
4. [Command Reference](#command-reference)
5. [Configuration](#configuration)
6. [Examples](#examples)

## Installation

### From Source

```bash
git clone https://github.com/dhineshp565/unified_pcr.git
cd unified_pcr
pip install -e .
```

### Requirements
- Python 3.7 or higher
- No additional dependencies for basic functionality

## Quick Start

Run the complete pipeline with a single command:

```bash
python -m src.ont_targseq.cli run \
    -i reads.fastq \
    -r targets.fasta \
    -o results/
```

This will:
1. Filter reads by quality and length
2. Align reads to target regions
3. Call variants (SNPs, indels)
4. Analyze coverage statistics

## Pipeline Steps

### 1. Preprocessing

Quality filter and trim adapter sequences from ONT reads.

```bash
python -m src.ont_targseq.cli preprocess \
    -i raw_reads.fastq \
    -o filtered_reads.fastq \
    --min-quality 7.0 \
    --min-length 100 \
    --max-length 10000 \
    --trim-adapters
```

**Options:**
- `--min-quality`: Minimum average Phred quality score (default: 7.0)
- `--min-length`: Minimum read length in bp (default: 100)
- `--max-length`: Maximum read length in bp (default: 10000)
- `--trim-adapters`: Enable adapter trimming

### 2. Alignment

Align filtered reads to target amplicon sequences.

```bash
python -m src.ont_targseq.cli align \
    -i filtered_reads.fastq \
    -r targets.fasta \
    -o alignments.sam \
    -t 4
```

**Options:**
- `-t, --threads`: Number of threads for alignment (default: 4)

### 3. Variant Calling

Identify SNPs, insertions, and deletions from aligned reads.

```bash
python -m src.ont_targseq.cli variants \
    -i alignments.sam \
    -r targets.fasta \
    -o variants.vcf \
    --min-coverage 10 \
    --min-frequency 0.2
```

**Options:**
- `--min-coverage`: Minimum read depth to call variants (default: 10)
- `--min-frequency`: Minimum allele frequency (0.0-1.0, default: 0.2)

### 4. Coverage Analysis

Analyze read coverage depth across target regions.

```bash
python -m src.ont_targseq.cli coverage \
    -i alignments.sam \
    -r targets.fasta \
    -o coverage.bed \
    --min-coverage 10
```

**Options:**
- `--min-coverage`: Threshold for coverage reporting (default: 10)

## Command Reference

### `run` - Full Pipeline

```bash
python -m src.ont_targseq.cli run \
    -i INPUT_FASTQ \
    -r REFERENCE_FASTA \
    -o OUTPUT_DIR \
    [-c CONFIG_FILE] \
    [--skip-preprocessing]
```

### `config` - Generate Config File

```bash
python -m src.ont_targseq.cli config -o config.json
```

This creates a default configuration file that you can edit and use with `run -c config.json`.

## Configuration

Configuration files use JSON format. Generate a template:

```bash
python -m src.ont_targseq.cli config -o my_config.json
```

Example configuration:

```json
{
  "preprocessing": {
    "min_quality": 7.0,
    "min_length": 100,
    "max_length": 10000,
    "trim_adapters": true,
    "adapter_sequences": []
  },
  "alignment": {
    "preset": "map-ont",
    "threads": 4,
    "min_mapping_quality": 20
  },
  "variant_calling": {
    "min_coverage": 10,
    "min_variant_frequency": 0.2,
    "min_base_quality": 7
  },
  "coverage_analysis": {
    "min_coverage": 10,
    "report_format": "bed"
  },
  "output": {
    "output_dir": "./ont_targseq_output",
    "keep_intermediates": true
  }
}
```

## Examples

### Example 1: Basic Analysis

```bash
# Run full pipeline with default settings
python -m src.ont_targseq.cli run \
    -i sample1.fastq \
    -r amplicons.fasta \
    -o sample1_results/
```

### Example 2: High Stringency Analysis

```bash
# Generate config with strict settings
python -m src.ont_targseq.cli config -o strict_config.json

# Edit config.json to increase min_quality to 10, min_coverage to 20, etc.

# Run with custom config
python -m src.ont_targseq.cli run \
    -i sample1.fastq \
    -r amplicons.fasta \
    -o sample1_strict/ \
    -c strict_config.json
```

### Example 3: Pre-filtered Data

```bash
# Skip preprocessing if reads are already filtered
python -m src.ont_targseq.cli run \
    -i pre_filtered.fastq \
    -r targets.fasta \
    -o results/ \
    --skip-preprocessing
```

### Example 4: Individual Steps

```bash
# Run each step separately for more control
python -m src.ont_targseq.cli preprocess -i raw.fastq -o filtered.fastq
python -m src.ont_targseq.cli align -i filtered.fastq -r targets.fasta -o align.sam
python -m src.ont_targseq.cli variants -i align.sam -r targets.fasta -o vars.vcf
python -m src.ont_targseq.cli coverage -i align.sam -r targets.fasta -o cov.bed
```

## Input File Formats

### FASTQ Format

Standard FASTQ format, can be plain text or gzipped (.fastq.gz):

```
@read_id
ACGTACGTACGT...
+
IIIIIIIIIIII...
```

### Reference FASTA Format

Target amplicon sequences in FASTA format:

```fasta
>BRCA1_exon2
ATCGATCGATCGATCG...
>BRCA2_exon5
GCTAGCTAGCTAGCTA...
```

## Output Files

| File | Description | Format |
|------|-------------|--------|
| filtered.fastq | Quality-filtered reads | FASTQ |
| alignments.sam | Aligned reads | SAM |
| variants.vcf | Called variants | VCF |
| coverage.bed | Coverage depth | BED |

## Troubleshooting

### Low Coverage

If you see low coverage warnings:
- Check input read quality
- Verify target sequences match your amplicons
- Reduce `--min-quality` threshold
- Check read length distribution

### No Variants Called

If no variants are detected:
- Lower `--min-frequency` threshold
- Reduce `--min-coverage` requirement
- Check alignment quality
- Verify reference sequences

### Memory Issues

For large datasets:
- Process samples individually
- Reduce thread count
- Split input FASTQ into smaller chunks

## Python API

For programmatic access:

```python
from src.ont_targseq import (
    FASTQPreprocessor,
    TargetAligner, 
    VariantCaller,
    CoverageAnalyzer
)

# Preprocess
preprocessor = FASTQPreprocessor(min_quality=7.0)
stats = preprocessor.filter_reads('input.fastq', 'filtered.fastq')

# Align
aligner = TargetAligner('reference.fasta')
aligner.align_reads('filtered.fastq', 'aligned.sam')

# Call variants
caller = VariantCaller(min_coverage=10)
caller.call_variants('aligned.sam', 'reference.fasta', 'variants.vcf')

# Analyze coverage
analyzer = CoverageAnalyzer(min_coverage=10)
coverage = analyzer.analyze_coverage('aligned.sam', 'reference.fasta')
```

## Performance Tips

1. **Use appropriate thread count**: `-t 8` for 8-core systems
2. **Pre-filter reads**: Remove obvious low-quality reads first
3. **Target specific regions**: Use focused reference FASTAs
4. **Batch processing**: Process multiple samples in parallel

## Getting Help

```bash
# General help
python -m src.ont_targseq.cli --help

# Command-specific help
python -m src.ont_targseq.cli preprocess --help
python -m src.ont_targseq.cli align --help
python -m src.ont_targseq.cli variants --help
python -m src.ont_targseq.cli coverage --help
python -m src.ont_targseq.cli run --help
```
