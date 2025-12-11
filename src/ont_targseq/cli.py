"""
Command-line interface for ONT target sequencing pipeline
"""

import argparse
import os
import sys
from typing import Optional

from .preprocessor import FASTQPreprocessor
from .aligner import TargetAligner
from .variant_caller import VariantCaller
from .coverage_analyzer import CoverageAnalyzer
from .config import Config


class ONTTargSeqCLI:
    """Command-line interface for the pipeline"""
    
    def __init__(self):
        self.parser = self._create_parser()
    
    def _create_parser(self) -> argparse.ArgumentParser:
        """Create argument parser"""
        parser = argparse.ArgumentParser(
            description='ONT Target Sequencing Pipeline - Unified PCR amplicon analysis',
            formatter_class=argparse.RawDescriptionHelpFormatter
        )
        
        subparsers = parser.add_subparsers(dest='command', help='Available commands')
        
        # Preprocess command
        preprocess = subparsers.add_parser('preprocess', help='Preprocess FASTQ files')
        preprocess.add_argument('-i', '--input', required=True, help='Input FASTQ file')
        preprocess.add_argument('-o', '--output', required=True, help='Output FASTQ file')
        preprocess.add_argument('--min-quality', type=float, default=7.0, help='Minimum quality score')
        preprocess.add_argument('--min-length', type=int, default=100, help='Minimum read length')
        preprocess.add_argument('--max-length', type=int, default=10000, help='Maximum read length')
        preprocess.add_argument('--trim-adapters', action='store_true', help='Trim adapter sequences')
        
        # Align command
        align = subparsers.add_parser('align', help='Align reads to targets')
        align.add_argument('-i', '--input', required=True, help='Input FASTQ file')
        align.add_argument('-r', '--reference', required=True, help='Reference FASTA file')
        align.add_argument('-o', '--output', required=True, help='Output SAM file')
        align.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
        
        # Call variants command
        variants = subparsers.add_parser('variants', help='Call variants from alignments')
        variants.add_argument('-i', '--input', required=True, help='Input SAM file')
        variants.add_argument('-r', '--reference', required=True, help='Reference FASTA file')
        variants.add_argument('-o', '--output', required=True, help='Output VCF file')
        variants.add_argument('--min-coverage', type=int, default=10, help='Minimum coverage')
        variants.add_argument('--min-frequency', type=float, default=0.2, help='Minimum variant frequency')
        
        # Coverage command
        coverage = subparsers.add_parser('coverage', help='Analyze coverage statistics')
        coverage.add_argument('-i', '--input', required=True, help='Input SAM file')
        coverage.add_argument('-r', '--reference', required=True, help='Reference FASTA file')
        coverage.add_argument('-o', '--output', help='Output BED file (optional)')
        coverage.add_argument('--min-coverage', type=int, default=10, help='Minimum coverage threshold')
        
        # Full pipeline command
        pipeline = subparsers.add_parser('run', help='Run full pipeline')
        pipeline.add_argument('-i', '--input', required=True, help='Input FASTQ file')
        pipeline.add_argument('-r', '--reference', required=True, help='Reference FASTA file')
        pipeline.add_argument('-o', '--output-dir', required=True, help='Output directory')
        pipeline.add_argument('-c', '--config', help='Configuration file (JSON)')
        pipeline.add_argument('--skip-preprocessing', action='store_true', help='Skip preprocessing step')
        
        # Config command
        config_cmd = subparsers.add_parser('config', help='Generate default config file')
        config_cmd.add_argument('-o', '--output', required=True, help='Output config file')
        
        return parser
    
    def run_preprocess(self, args):
        """Run preprocessing step"""
        print(f"Preprocessing reads from {args.input}...")
        
        preprocessor = FASTQPreprocessor(
            min_quality=args.min_quality,
            min_length=args.min_length,
            max_length=args.max_length
        )
        
        stats = preprocessor.filter_reads(args.input, args.output)
        
        print(f"Total reads: {stats['total_reads']}")
        print(f"Passed reads: {stats['passed_reads']}")
        print(f"Failed quality: {stats['failed_quality']}")
        print(f"Failed length: {stats['failed_length']}")
        
        if args.trim_adapters:
            print("\nTrimming adapters...")
            trimmed_output = args.output.replace('.fastq', '_trimmed.fastq')
            trim_stats = preprocessor.trim_adapters(args.output, trimmed_output)
            print(f"Trimmed reads: {trim_stats['trimmed_reads']}")
            print(f"Bases trimmed: {trim_stats['bases_trimmed']}")
    
    def run_align(self, args):
        """Run alignment step"""
        print(f"Aligning reads to {args.reference}...")
        
        aligner = TargetAligner(args.reference)
        targets = aligner.get_target_regions()
        
        print(f"Loaded {len(targets)} target regions")
        for name, (seq, length) in targets.items():
            print(f"  {name}: {length} bp")
        
        stats = aligner.align_reads(args.input, args.output, threads=args.threads)
        print(f"\nAlignment complete. Output saved to {args.output}")
    
    def run_variants(self, args):
        """Run variant calling step"""
        print(f"Calling variants from {args.input}...")
        
        caller = VariantCaller(
            min_coverage=args.min_coverage,
            min_variant_frequency=args.min_frequency
        )
        
        stats = caller.call_variants(args.input, args.reference, args.output)
        
        print(f"Total positions analyzed: {stats['total_positions']}")
        print(f"Variant positions: {stats['variant_positions']}")
        print(f"SNPs: {stats['snps']}")
        print(f"Insertions: {stats['insertions']}")
        print(f"Deletions: {stats['deletions']}")
        print(f"\nVariants written to {args.output}")
    
    def run_coverage(self, args):
        """Run coverage analysis step"""
        print(f"Analyzing coverage from {args.input}...")
        
        analyzer = CoverageAnalyzer(min_coverage=args.min_coverage)
        stats = analyzer.analyze_coverage(args.input, args.reference, args.output)
        
        print("\n=== Overall Coverage Statistics ===")
        overall = stats['overall']
        print(f"Mean coverage: {overall['mean_coverage']:.2f}x")
        print(f"Median coverage: {overall['median_coverage']}x")
        print(f"Min coverage: {overall['min_coverage']}x")
        print(f"Max coverage: {overall['max_coverage']}x")
        print(f"Bases above {args.min_coverage}x: {overall['percent_above_threshold']:.2f}%")
        
        print("\n=== Per-Target Coverage ===")
        for target, target_stats in stats['targets'].items():
            print(f"\n{target}:")
            print(f"  Length: {target_stats['length']} bp")
            print(f"  Mean: {target_stats['mean_coverage']:.2f}x")
            print(f"  Median: {target_stats['median_coverage']}x")
            print(f"  Coverage >= {args.min_coverage}x: {target_stats['percent_above_threshold']:.2f}%")
        
        if args.output:
            print(f"\nCoverage data written to {args.output}")
    
    def run_pipeline(self, args):
        """Run full pipeline"""
        print("=== ONT Target Sequencing Pipeline ===\n")
        
        # Load configuration
        config = Config(args.config) if args.config else Config()
        
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Step 1: Preprocessing
        if not args.skip_preprocessing:
            print("Step 1: Preprocessing...")
            filtered_fastq = os.path.join(args.output_dir, 'filtered.fastq')
            preprocess_args = argparse.Namespace(
                input=args.input,
                output=filtered_fastq,
                min_quality=config.get('preprocessing', 'min_quality'),
                min_length=config.get('preprocessing', 'min_length'),
                max_length=config.get('preprocessing', 'max_length'),
                trim_adapters=config.get('preprocessing', 'trim_adapters')
            )
            self.run_preprocess(preprocess_args)
            input_fastq = filtered_fastq
        else:
            input_fastq = args.input
        
        # Step 2: Alignment
        print("\nStep 2: Alignment...")
        sam_file = os.path.join(args.output_dir, 'alignments.sam')
        align_args = argparse.Namespace(
            input=input_fastq,
            reference=args.reference,
            output=sam_file,
            threads=config.get('alignment', 'threads')
        )
        self.run_align(align_args)
        
        # Step 3: Variant calling
        print("\nStep 3: Variant calling...")
        vcf_file = os.path.join(args.output_dir, 'variants.vcf')
        variant_args = argparse.Namespace(
            input=sam_file,
            reference=args.reference,
            output=vcf_file,
            min_coverage=config.get('variant_calling', 'min_coverage'),
            min_frequency=config.get('variant_calling', 'min_variant_frequency')
        )
        self.run_variants(variant_args)
        
        # Step 4: Coverage analysis
        print("\nStep 4: Coverage analysis...")
        bed_file = os.path.join(args.output_dir, 'coverage.bed')
        coverage_args = argparse.Namespace(
            input=sam_file,
            reference=args.reference,
            output=bed_file,
            min_coverage=config.get('coverage_analysis', 'min_coverage')
        )
        self.run_coverage(coverage_args)
        
        print(f"\n=== Pipeline Complete ===")
        print(f"Results saved to: {args.output_dir}")
    
    def run_config(self, args):
        """Generate default configuration file"""
        config = Config()
        config.create_default_config(args.output)
        print(f"Default configuration written to {args.output}")
    
    def run(self, argv=None):
        """Run the CLI"""
        args = self.parser.parse_args(argv)
        
        if not args.command:
            self.parser.print_help()
            return 1
        
        try:
            if args.command == 'preprocess':
                self.run_preprocess(args)
            elif args.command == 'align':
                self.run_align(args)
            elif args.command == 'variants':
                self.run_variants(args)
            elif args.command == 'coverage':
                self.run_coverage(args)
            elif args.command == 'run':
                self.run_pipeline(args)
            elif args.command == 'config':
                self.run_config(args)
            else:
                self.parser.print_help()
                return 1
            
            return 0
        except Exception as e:
            print(f"Error: {str(e)}", file=sys.stderr)
            return 1


def main():
    """Main entry point"""
    cli = ONTTargSeqCLI()
    sys.exit(cli.run())


if __name__ == '__main__':
    main()
