# unified_pcr

Pipeline for targeted amplicon sequencing using Oxford Nanopore Technologies (ONT)

## Overview

This repository contains a comprehensive pipeline for analyzing ONT targeted amplicon sequencing data. The pipeline includes modules for:

- FASTQ preprocessing and quality filtering
- Target region alignment
- Variant calling (SNPs, insertions, deletions)
- Coverage analysis

## Quick Start

```bash
# Install the package
pip install -e .

# Run the full pipeline
python -m src.ont_targseq.cli run \
    -i input_reads.fastq \
    -r targets.fasta \
    -o output_directory
```

## Documentation

See [docs/README.md](docs/README.md) for detailed documentation.

## Examples

Check out the [examples/](examples/) directory for example workflows and usage patterns.

## Testing

Run tests with:
```bash
python -m pytest tests/
```

## License

See LICENSE file for details.
