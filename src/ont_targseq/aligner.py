"""
Target Aligner for ONT reads
Aligns reads to target regions using minimap2-style alignment
"""

import os
from typing import Dict, List, Optional, Tuple


class TargetAligner:
    """Aligns ONT reads to target reference sequences"""
    
    def __init__(self, reference_fasta: str, preset: str = 'map-ont'):
        """
        Initialize the target aligner
        
        Args:
            reference_fasta: Path to reference FASTA file
            preset: Alignment preset (map-ont for Oxford Nanopore)
        """
        self.reference_fasta = reference_fasta
        self.preset = preset
        self.targets = {}
        self._load_reference()
    
    def _load_reference(self):
        """Load reference sequences from FASTA file"""
        if not os.path.exists(self.reference_fasta):
            return
            
        current_seq_name = None
        current_seq = []
        
        with open(self.reference_fasta, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq_name:
                        self.targets[current_seq_name] = ''.join(current_seq)
                    current_seq_name = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            
            if current_seq_name:
                self.targets[current_seq_name] = ''.join(current_seq)
    
    def align_reads(self, input_fastq: str, output_sam: str, 
                    threads: int = 4) -> Dict[str, int]:
        """
        Align reads to target regions
        
        Args:
            input_fastq: Path to input FASTQ file
            output_sam: Path to output SAM file
            threads: Number of threads for alignment
            
        Returns:
            Dictionary with alignment statistics
        """
        stats = {
            'total_reads': 0,
            'aligned_reads': 0,
            'unaligned_reads': 0,
            'primary_alignments': 0
        }
        
        # Simplified alignment - in practice would use minimap2 or similar
        # This is a placeholder implementation
        with open(output_sam, 'w') as sam_out:
            # Write SAM header
            sam_out.write("@HD\tVN:1.6\tSO:unsorted\n")
            for target_name, target_seq in self.targets.items():
                sam_out.write(f"@SQ\tSN:{target_name}\tLN:{len(target_seq)}\n")
            
            # In a real implementation, this would call minimap2 or similar aligner
            # For now, just write a placeholder comment
            sam_out.write(f"@PG\tID:ont_targseq\tPN:ont_targseq\tVN:1.0.0\n")
        
        return stats
    
    def get_target_regions(self) -> Dict[str, Tuple[str, int]]:
        """
        Get information about target regions
        
        Returns:
            Dictionary mapping target names to (sequence, length)
        """
        return {name: (seq, len(seq)) for name, seq in self.targets.items()}
    
    def simple_align(self, query_seq: str, target_name: str, 
                     min_match_length: int = 20) -> Optional[Dict]:
        """
        Simple alignment of a query sequence to a target
        
        Args:
            query_seq: Query sequence to align
            target_name: Name of target sequence
            min_match_length: Minimum match length to report
            
        Returns:
            Alignment information dictionary or None
        """
        if target_name not in self.targets:
            return None
        
        target_seq = self.targets[target_name]
        
        # Simple substring search (in practice use proper alignment)
        best_match = None
        for i in range(len(target_seq) - min_match_length + 1):
            for j in range(min_match_length, min(len(query_seq), len(target_seq) - i) + 1):
                target_substring = target_seq[i:i+j]
                if target_substring in query_seq:
                    if best_match is None or j > best_match['length']:
                        best_match = {
                            'target_start': i,
                            'target_end': i + j,
                            'length': j,
                            'query_match': target_substring
                        }
        
        return best_match
