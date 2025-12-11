"""
Variant Caller for ONT target sequencing data
Identifies variants in aligned reads
"""

import os
from typing import Dict, List, Optional, Tuple
from collections import defaultdict


class VariantCaller:
    """Identifies variants from aligned ONT reads"""
    
    def __init__(self, min_coverage: int = 10, min_variant_frequency: float = 0.2):
        """
        Initialize the variant caller
        
        Args:
            min_coverage: Minimum coverage required to call a variant
            min_variant_frequency: Minimum allele frequency to report variant
        """
        self.min_coverage = min_coverage
        self.min_variant_frequency = min_variant_frequency
    
    def call_variants(self, sam_file: str, reference_fasta: str, 
                      output_vcf: str) -> Dict[str, int]:
        """
        Call variants from SAM alignment file
        
        Args:
            sam_file: Path to SAM alignment file
            reference_fasta: Path to reference FASTA file
            output_vcf: Path to output VCF file
            
        Returns:
            Dictionary with variant calling statistics
        """
        stats = {
            'total_positions': 0,
            'variant_positions': 0,
            'snps': 0,
            'insertions': 0,
            'deletions': 0
        }
        
        # Load reference
        reference = self._load_reference(reference_fasta)
        
        # Parse alignments and build pileup
        pileup = self._build_pileup(sam_file)
        
        # Call variants
        variants = []
        for chrom in pileup:
            for pos in sorted(pileup[chrom].keys()):
                stats['total_positions'] += 1
                bases = pileup[chrom][pos]
                coverage = sum(bases.values())
                
                if coverage < self.min_coverage:
                    continue
                
                # Get reference base
                ref_seq = reference.get(chrom, '')
                ref_base = ref_seq[pos] if pos < len(ref_seq) else 'N'
                
                # Find variant alleles
                for alt_base, count in bases.items():
                    if alt_base != ref_base:
                        freq = count / coverage
                        if freq >= self.min_variant_frequency:
                            variants.append({
                                'chrom': chrom,
                                'pos': pos + 1,  # VCF is 1-based
                                'ref': ref_base,
                                'alt': alt_base,
                                'qual': int(freq * 100),
                                'coverage': coverage,
                                'alt_count': count,
                                'frequency': freq
                            })
                            stats['variant_positions'] += 1
                            
                            # Classify variant type
                            if len(ref_base) == 1 and len(alt_base) == 1:
                                stats['snps'] += 1
                            elif len(alt_base) > len(ref_base):
                                stats['insertions'] += 1
                            else:
                                stats['deletions'] += 1
        
        # Write VCF file
        self._write_vcf(variants, output_vcf, reference)
        
        return stats
    
    def _load_reference(self, fasta_file: str) -> Dict[str, str]:
        """Load reference sequences from FASTA"""
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
    
    def _build_pileup(self, sam_file: str) -> Dict[str, Dict[int, Dict[str, int]]]:
        """Build pileup from SAM file"""
        pileup = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        
        if not os.path.exists(sam_file):
            return pileup
        
        with open(sam_file, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 11:
                    continue
                
                chrom = fields[2]
                pos = int(fields[3]) - 1  # Convert to 0-based
                seq = fields[9]
                
                # Simple pileup (ignores CIGAR in this simplified version)
                for i, base in enumerate(seq):
                    pileup[chrom][pos + i][base] += 1
        
        return pileup
    
    def _write_vcf(self, variants: List[Dict], output_vcf: str, 
                   reference: Dict[str, str]):
        """Write variants to VCF file"""
        with open(output_vcf, 'w') as f:
            # Write VCF header
            f.write("##fileformat=VCFv4.2\n")
            f.write("##source=ont_targseq_v1.0.0\n")
            
            # Write reference contigs
            for chrom, seq in reference.items():
                f.write(f"##contig=<ID={chrom},length={len(seq)}>\n")
            
            # Write INFO and FORMAT headers
            f.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
            f.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
            f.write("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele Count\">\n")
            
            # Write column headers
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            
            # Write variants
            for var in variants:
                info = f"DP={var['coverage']};AF={var['frequency']:.3f};AC={var['alt_count']}"
                f.write(f"{var['chrom']}\t{var['pos']}\t.\t{var['ref']}\t{var['alt']}\t"
                       f"{var['qual']}\tPASS\t{info}\n")
