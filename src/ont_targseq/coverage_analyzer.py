"""
Coverage Analyzer for ONT target sequencing
Analyzes coverage depth and uniformity across target regions
"""

import os
from typing import Dict, List, Optional, Tuple
from collections import defaultdict


class CoverageAnalyzer:
    """Analyzes coverage statistics for target regions"""
    
    def __init__(self, min_coverage: int = 10):
        """
        Initialize the coverage analyzer
        
        Args:
            min_coverage: Minimum coverage threshold for reporting
        """
        self.min_coverage = min_coverage
    
    def analyze_coverage(self, sam_file: str, reference_fasta: str,
                        output_bed: Optional[str] = None) -> Dict[str, any]:
        """
        Analyze coverage across target regions
        
        Args:
            sam_file: Path to SAM alignment file
            reference_fasta: Path to reference FASTA file
            output_bed: Optional path to output BED file with coverage
            
        Returns:
            Dictionary with coverage statistics
        """
        # Load reference to get target lengths
        reference = self._load_reference(reference_fasta)
        
        # Build coverage array
        coverage = self._calculate_coverage(sam_file, reference)
        
        # Calculate statistics
        stats = self._calculate_stats(coverage, reference)
        
        # Write BED file if requested
        if output_bed:
            self._write_bed(coverage, output_bed)
        
        return stats
    
    def _load_reference(self, fasta_file: str) -> Dict[str, str]:
        """Load reference sequences"""
        reference = {}
        if not os.path.exists(fasta_file):
            return reference
            
        current_chrom = None
        current_seq = []
        
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_chrom:
                        reference[current_chrom] = ''.join(current_seq)
                    current_chrom = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            
            if current_chrom:
                reference[current_chrom] = ''.join(current_seq)
        
        return reference
    
    def _calculate_coverage(self, sam_file: str, 
                           reference: Dict[str, str]) -> Dict[str, List[int]]:
        """Calculate per-base coverage from SAM file"""
        coverage = {}
        
        # Initialize coverage arrays
        for chrom, seq in reference.items():
            coverage[chrom] = [0] * len(seq)
        
        if not os.path.exists(sam_file):
            return coverage
        
        # Parse SAM file and count coverage
        with open(sam_file, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 11:
                    continue
                
                chrom = fields[2]
                if chrom not in coverage:
                    continue
                
                pos = int(fields[3]) - 1  # Convert to 0-based
                seq = fields[9]
                
                # Simple coverage calculation (ignores CIGAR)
                for i in range(len(seq)):
                    if pos + i < len(coverage[chrom]):
                        coverage[chrom][pos + i] += 1
        
        return coverage
    
    def _calculate_stats(self, coverage: Dict[str, List[int]], 
                        reference: Dict[str, str]) -> Dict[str, any]:
        """Calculate coverage statistics"""
        stats = {
            'targets': {},
            'overall': {
                'mean_coverage': 0.0,
                'median_coverage': 0.0,
                'min_coverage': 0,
                'max_coverage': 0,
                'bases_above_threshold': 0,
                'total_bases': 0,
                'percent_above_threshold': 0.0
            }
        }
        
        all_coverages = []
        total_coverage = 0
        bases_above_threshold = 0
        total_bases = 0
        
        for chrom in coverage:
            chrom_cov = coverage[chrom]
            if not chrom_cov:
                continue
            
            # Calculate per-target statistics
            mean_cov = sum(chrom_cov) / len(chrom_cov) if chrom_cov else 0
            sorted_cov = sorted(chrom_cov)
            median_cov = sorted_cov[len(sorted_cov) // 2] if sorted_cov else 0
            min_cov = min(chrom_cov) if chrom_cov else 0
            max_cov = max(chrom_cov) if chrom_cov else 0
            above_threshold = sum(1 for c in chrom_cov if c >= self.min_coverage)
            
            stats['targets'][chrom] = {
                'length': len(chrom_cov),
                'mean_coverage': mean_cov,
                'median_coverage': median_cov,
                'min_coverage': min_cov,
                'max_coverage': max_cov,
                'bases_above_threshold': above_threshold,
                'percent_above_threshold': (above_threshold / len(chrom_cov) * 100) if chrom_cov else 0
            }
            
            # Aggregate for overall stats
            all_coverages.extend(chrom_cov)
            total_coverage += sum(chrom_cov)
            total_bases += len(chrom_cov)
            bases_above_threshold += above_threshold
        
        # Calculate overall statistics
        if all_coverages:
            stats['overall']['mean_coverage'] = total_coverage / total_bases if total_bases > 0 else 0
            sorted_all = sorted(all_coverages)
            stats['overall']['median_coverage'] = sorted_all[len(sorted_all) // 2]
            stats['overall']['min_coverage'] = min(all_coverages)
            stats['overall']['max_coverage'] = max(all_coverages)
            stats['overall']['bases_above_threshold'] = bases_above_threshold
            stats['overall']['total_bases'] = total_bases
            stats['overall']['percent_above_threshold'] = (bases_above_threshold / total_bases * 100) if total_bases > 0 else 0
        
        return stats
    
    def _write_bed(self, coverage: Dict[str, List[int]], output_bed: str):
        """Write coverage to BED file"""
        with open(output_bed, 'w') as f:
            for chrom in sorted(coverage.keys()):
                for pos, cov in enumerate(coverage[chrom]):
                    f.write(f"{chrom}\t{pos}\t{pos+1}\t{cov}\n")
    
    def get_low_coverage_regions(self, coverage_stats: Dict, 
                                 threshold: Optional[int] = None) -> List[Dict]:
        """
        Identify regions with low coverage
        
        Args:
            coverage_stats: Output from analyze_coverage
            threshold: Coverage threshold (uses min_coverage if not specified)
            
        Returns:
            List of low coverage regions
        """
        if threshold is None:
            threshold = self.min_coverage
        
        low_coverage_regions = []
        
        for target, stats in coverage_stats.get('targets', {}).items():
            if stats['mean_coverage'] < threshold:
                low_coverage_regions.append({
                    'target': target,
                    'mean_coverage': stats['mean_coverage'],
                    'length': stats['length']
                })
        
        return low_coverage_regions
