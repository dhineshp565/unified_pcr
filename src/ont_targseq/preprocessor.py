"""
FASTQ Preprocessor for ONT data
Handles quality filtering and adapter trimming for Oxford Nanopore reads
"""

import os
from typing import Dict, List, Optional, Tuple
import gzip


class FASTQPreprocessor:
    """Preprocessor for ONT FASTQ files with quality filtering"""
    
    def __init__(self, min_quality: float = 7.0, min_length: int = 100, max_length: int = 10000):
        """
        Initialize the FASTQ preprocessor
        
        Args:
            min_quality: Minimum average quality score (Phred)
            min_length: Minimum read length to keep
            max_length: Maximum read length to keep
        """
        self.min_quality = min_quality
        self.min_length = min_length
        self.max_length = max_length
        
    def _calculate_average_quality(self, quality_string: str) -> float:
        """Calculate average Phred quality score from quality string"""
        if not quality_string:
            return 0.0
        return sum(ord(char) - 33 for char in quality_string) / len(quality_string)
    
    def _open_fastq(self, filepath: str):
        """Open FASTQ file (handles gzipped files)"""
        if filepath.endswith('.gz'):
            return gzip.open(filepath, 'rt')
        return open(filepath, 'r')
    
    def filter_reads(self, input_fastq: str, output_fastq: str) -> Dict[str, int]:
        """
        Filter reads based on quality and length criteria
        
        Args:
            input_fastq: Path to input FASTQ file
            output_fastq: Path to output filtered FASTQ file
            
        Returns:
            Dictionary with filtering statistics
        """
        stats = {
            'total_reads': 0,
            'passed_reads': 0,
            'failed_quality': 0,
            'failed_length': 0
        }
        
        with self._open_fastq(input_fastq) as infile, open(output_fastq, 'w') as outfile:
            while True:
                # Read four lines (one FASTQ record)
                header = infile.readline().strip()
                if not header:
                    break
                    
                sequence = infile.readline().strip()
                plus = infile.readline().strip()
                quality = infile.readline().strip()
                
                stats['total_reads'] += 1
                
                # Check length
                seq_length = len(sequence)
                if seq_length < self.min_length or seq_length > self.max_length:
                    stats['failed_length'] += 1
                    continue
                
                # Check quality
                avg_quality = self._calculate_average_quality(quality)
                if avg_quality < self.min_quality:
                    stats['failed_quality'] += 1
                    continue
                
                # Write passing read
                outfile.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")
                stats['passed_reads'] += 1
        
        return stats
    
    def trim_adapters(self, input_fastq: str, output_fastq: str, 
                      adapter_sequences: Optional[List[str]] = None) -> Dict[str, int]:
        """
        Trim adapter sequences from reads
        
        Args:
            input_fastq: Path to input FASTQ file
            output_fastq: Path to output trimmed FASTQ file
            adapter_sequences: List of adapter sequences to trim
            
        Returns:
            Dictionary with trimming statistics
        """
        if adapter_sequences is None:
            # Default ONT adapters
            adapter_sequences = [
                'AATGTACTTCGTTCAGTTACGTATTGCT',  # Example ONT adapter
            ]
        
        stats = {
            'total_reads': 0,
            'trimmed_reads': 0,
            'bases_trimmed': 0
        }
        
        with self._open_fastq(input_fastq) as infile, open(output_fastq, 'w') as outfile:
            while True:
                header = infile.readline().strip()
                if not header:
                    break
                    
                sequence = infile.readline().strip()
                plus = infile.readline().strip()
                quality = infile.readline().strip()
                
                stats['total_reads'] += 1
                original_length = len(sequence)
                
                # Simple adapter trimming (look for adapter at end)
                trimmed = False
                for adapter in adapter_sequences:
                    if adapter in sequence:
                        pos = sequence.find(adapter)
                        if pos > 0:
                            sequence = sequence[:pos]
                            quality = quality[:pos]
                            trimmed = True
                            break
                
                if trimmed:
                    stats['trimmed_reads'] += 1
                    stats['bases_trimmed'] += original_length - len(sequence)
                
                outfile.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")
        
        return stats
