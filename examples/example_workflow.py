"""
Example workflow for ONT target sequencing pipeline
"""

import os
import sys

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.ont_targseq import FASTQPreprocessor, TargetAligner, VariantCaller, CoverageAnalyzer
from src.ont_targseq.config import Config


def create_example_data(output_dir='example_data'):
    """Create example input files for demonstration"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Create example reference FASTA
    reference_fasta = os.path.join(output_dir, 'targets.fasta')
    with open(reference_fasta, 'w') as f:
        f.write('>GENE1_exon1\n')
        f.write('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n')
        f.write('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n')
        f.write('>GENE2_exon2\n')
        f.write('GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n')
        f.write('GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n')
    
    # Create example FASTQ
    input_fastq = os.path.join(output_dir, 'reads.fastq')
    with open(input_fastq, 'w') as f:
        # High quality read
        f.write('@read1\n')
        f.write('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n')
        f.write('+\n')
        f.write('IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n')
        
        # Low quality read
        f.write('@read2\n')
        f.write('GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n')
        f.write('+\n')
        f.write('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
        
        # Short read
        f.write('@read3\n')
        f.write('ATCGATCG\n')
        f.write('+\n')
        f.write('IIIIIIII\n')
        
        # Good read
        f.write('@read4\n')
        f.write('GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n')
        f.write('+\n')
        f.write('IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n')
    
    return reference_fasta, input_fastq


def run_example_workflow():
    """Run an example workflow"""
    print("=== ONT Target Sequencing Example Workflow ===\n")
    
    # Create example data
    print("Creating example data...")
    reference_fasta, input_fastq = create_example_data()
    print(f"Reference: {reference_fasta}")
    print(f"Input FASTQ: {input_fastq}")
    
    # Create output directory
    output_dir = 'example_output'
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Preprocessing
    print("\n--- Step 1: Preprocessing ---")
    preprocessor = FASTQPreprocessor(min_quality=7.0, min_length=50, max_length=10000)
    filtered_fastq = os.path.join(output_dir, 'filtered.fastq')
    stats = preprocessor.filter_reads(input_fastq, filtered_fastq)
    
    print(f"Total reads: {stats['total_reads']}")
    print(f"Passed reads: {stats['passed_reads']}")
    print(f"Failed quality: {stats['failed_quality']}")
    print(f"Failed length: {stats['failed_length']}")
    print(f"Output: {filtered_fastq}")
    
    # Step 2: Alignment
    print("\n--- Step 2: Alignment ---")
    aligner = TargetAligner(reference_fasta, preset='map-ont')
    targets = aligner.get_target_regions()
    print(f"Loaded {len(targets)} target regions:")
    for name, (seq, length) in targets.items():
        print(f"  {name}: {length} bp")
    
    sam_file = os.path.join(output_dir, 'alignments.sam')
    align_stats = aligner.align_reads(filtered_fastq, sam_file, threads=4)
    print(f"Output: {sam_file}")
    
    # Step 3: Variant calling
    print("\n--- Step 3: Variant Calling ---")
    caller = VariantCaller(min_coverage=5, min_variant_frequency=0.1)
    vcf_file = os.path.join(output_dir, 'variants.vcf')
    variant_stats = caller.call_variants(sam_file, reference_fasta, vcf_file)
    
    print(f"Total positions analyzed: {variant_stats['total_positions']}")
    print(f"Variant positions: {variant_stats['variant_positions']}")
    print(f"SNPs: {variant_stats['snps']}")
    print(f"Insertions: {variant_stats['insertions']}")
    print(f"Deletions: {variant_stats['deletions']}")
    print(f"Output: {vcf_file}")
    
    # Step 4: Coverage analysis
    print("\n--- Step 4: Coverage Analysis ---")
    analyzer = CoverageAnalyzer(min_coverage=5)
    bed_file = os.path.join(output_dir, 'coverage.bed')
    coverage_stats = analyzer.analyze_coverage(sam_file, reference_fasta, bed_file)
    
    print("\nOverall Statistics:")
    overall = coverage_stats['overall']
    print(f"  Mean coverage: {overall['mean_coverage']:.2f}x")
    print(f"  Median coverage: {overall['median_coverage']}x")
    print(f"  Coverage >= 5x: {overall['percent_above_threshold']:.2f}%")
    
    print("\nPer-Target Statistics:")
    for target, target_stats in coverage_stats['targets'].items():
        print(f"\n  {target}:")
        print(f"    Length: {target_stats['length']} bp")
        print(f"    Mean coverage: {target_stats['mean_coverage']:.2f}x")
        print(f"    Coverage >= 5x: {target_stats['percent_above_threshold']:.2f}%")
    
    print(f"\nOutput: {bed_file}")
    
    # Summary
    print("\n=== Workflow Complete ===")
    print(f"All output files saved to: {output_dir}/")
    print("\nGenerated files:")
    for filename in os.listdir(output_dir):
        filepath = os.path.join(output_dir, filename)
        size = os.path.getsize(filepath)
        print(f"  {filename} ({size} bytes)")


def run_with_config():
    """Run workflow with configuration file"""
    print("\n=== Example: Using Configuration File ===\n")
    
    # Create and save configuration
    config = Config()
    config.set('preprocessing', 'min_quality', 10.0)
    config.set('variant_calling', 'min_coverage', 15)
    
    config_file = 'example_config.json'
    config.save_config(config_file)
    print(f"Configuration saved to: {config_file}")
    
    # Load and use configuration
    loaded_config = Config(config_file)
    print(f"Min quality: {loaded_config.get('preprocessing', 'min_quality')}")
    print(f"Min coverage: {loaded_config.get('variant_calling', 'min_coverage')}")


if __name__ == '__main__':
    # Run the example workflow
    run_example_workflow()
    
    # Show configuration example
    run_with_config()
