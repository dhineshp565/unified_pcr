# ONT Target Sequencing Pipeline

A unified pipeline for Oxford Nanopore Technologies (ONT) targeted amplicon sequencing analysis.

## Features

- **FASTQ Preprocessing**: Quality filtering and adapter trimming for ONT reads
- **Target Alignment**: Align reads to target amplicon regions
- **Variant Calling**: Identify SNPs, insertions, and deletions
- **Coverage Analysis**: Analyze coverage depth and uniformity across targets

## Installation

```bash
# Clone the repository
git clone https://github.com/dhineshp565/unified_pcr.git
cd unified_pcr

# Install dependencies
pip install -r requirements.txt

# Install the package
pip install -e .
```

## Quick Start

### Run Full Pipeline

```bash
python -m src.ont_targseq.cli run \
    -i input_reads.fastq \
    -r targets.fasta \
    -o output_directory
```

### Individual Steps

#### 1. Preprocess Reads

```bash
python -m src.ont_targseq.cli preprocess \
    -i raw_reads.fastq \
    -o filtered_reads.fastq \
    --min-quality 7.0 \
    --min-length 100 \
    --trim-adapters
```

#### 2. Align to Targets

```bash
python -m src.ont_targseq.cli align \
    -i filtered_reads.fastq \
    -r targets.fasta \
    -o alignments.sam \
    -t 4
```

#### 3. Call Variants

```bash
python -m src.ont_targseq.cli variants \
    -i alignments.sam \
    -r targets.fasta \
    -o variants.vcf \
    --min-coverage 10 \
    --min-frequency 0.2
```

#### 4. Analyze Coverage

```bash
python -m src.ont_targseq.cli coverage \
    -i alignments.sam \
    -r targets.fasta \
    -o coverage.bed \
    --min-coverage 10
```

## Configuration

Generate a default configuration file:

```bash
python -m src.ont_targseq.cli config -o config.json
```

Edit the configuration file to customize pipeline parameters, then run:

```bash
python -m src.ont_targseq.cli run \
    -i input_reads.fastq \
    -r targets.fasta \
    -o output_directory \
    -c config.json
```

### Configuration Options

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

## Python API

```python
from src.ont_targseq import FASTQPreprocessor, TargetAligner, VariantCaller, CoverageAnalyzer

# Preprocess reads
preprocessor = FASTQPreprocessor(min_quality=7.0, min_length=100)
stats = preprocessor.filter_reads('input.fastq', 'filtered.fastq')

# Align reads
aligner = TargetAligner('reference.fasta')
aligner.align_reads('filtered.fastq', 'alignments.sam')

# Call variants
caller = VariantCaller(min_coverage=10, min_variant_frequency=0.2)
caller.call_variants('alignments.sam', 'reference.fasta', 'variants.vcf')

# Analyze coverage
analyzer = CoverageAnalyzer(min_coverage=10)
coverage_stats = analyzer.analyze_coverage('alignments.sam', 'reference.fasta', 'coverage.bed')
```

## Input Format

### FASTQ File
Standard FASTQ format with ONT reads. Supports both plain and gzipped files (.fastq.gz).

### Reference FASTA File
FASTA file containing target amplicon sequences:

```
>target1
ACGTACGTACGTACGTACGTACGTACGTACGTACGT
>target2
TGCATGCATGCATGCATGCATGCATGCATGCA
```

## Output Files

- **filtered.fastq**: Quality-filtered reads
- **alignments.sam**: Aligned reads in SAM format
- **variants.vcf**: Called variants in VCF format
- **coverage.bed**: Coverage depth in BED format

## Testing

Run the test suite:

```bash
python -m pytest tests/
```

Or using unittest:

```bash
python -m unittest discover tests/
```

## Pipeline Workflow

```
Raw FASTQ
    ↓
[Preprocessing]
    ├─ Quality filtering
    └─ Adapter trimming
    ↓
Filtered FASTQ
    ↓
[Alignment]
    └─ Map to target regions
    ↓
SAM alignments
    ├─ [Variant Calling]
    │      └─ VCF output
    └─ [Coverage Analysis]
           └─ BED output
```

## Requirements

- Python 3.7+
- See requirements.txt for Python dependencies

## License

See LICENSE file for details.

## Support

For issues and questions, please open an issue on GitHub.
