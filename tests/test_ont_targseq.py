"""
Unit tests for ONT target sequencing pipeline
"""

import unittest
import os
import tempfile
import shutil
from src.ont_targseq.preprocessor import FASTQPreprocessor
from src.ont_targseq.aligner import TargetAligner
from src.ont_targseq.variant_caller import VariantCaller
from src.ont_targseq.coverage_analyzer import CoverageAnalyzer
from src.ont_targseq.config import Config


class TestFASTQPreprocessor(unittest.TestCase):
    """Test FASTQ preprocessing functionality"""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.test_fastq = os.path.join(self.temp_dir, 'test.fastq')
        
        # Create a simple test FASTQ file
        with open(self.test_fastq, 'w') as f:
            f.write('@read1\n')
            f.write('ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n')
            f.write('+\n')
            f.write('IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n')
            f.write('@read2\n')
            f.write('ACGT\n')  # Too short
            f.write('+\n')
            f.write('IIII\n')
            f.write('@read3\n')
            f.write('ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n')
            f.write('+\n')
            f.write('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')  # Low quality
    
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
    
    def test_filter_reads(self):
        """Test read filtering"""
        preprocessor = FASTQPreprocessor(min_quality=7.0, min_length=10)
        output_fastq = os.path.join(self.temp_dir, 'filtered.fastq')
        
        stats = preprocessor.filter_reads(self.test_fastq, output_fastq)
        
        self.assertEqual(stats['total_reads'], 3)
        self.assertEqual(stats['passed_reads'], 1)
        self.assertEqual(stats['failed_length'], 1)
        self.assertEqual(stats['failed_quality'], 1)
    
    def test_quality_calculation(self):
        """Test average quality calculation"""
        preprocessor = FASTQPreprocessor()
        
        # 'I' has Phred quality of ord('I') - 33 = 73 - 33 = 40
        quality = preprocessor._calculate_average_quality('IIII')
        self.assertEqual(quality, 40.0)
        
        # '!' has Phred quality of ord('!') - 33 = 33 - 33 = 0
        quality = preprocessor._calculate_average_quality('!!!!')
        self.assertEqual(quality, 0.0)


class TestTargetAligner(unittest.TestCase):
    """Test target alignment functionality"""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.test_reference = os.path.join(self.temp_dir, 'reference.fasta')
        
        # Create a simple reference FASTA
        with open(self.test_reference, 'w') as f:
            f.write('>target1\n')
            f.write('ACGTACGTACGTACGTACGTACGTACGTACGT\n')
            f.write('>target2\n')
            f.write('TGCATGCATGCATGCATGCATGCATGCATGCA\n')
    
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
    
    def test_load_reference(self):
        """Test reference loading"""
        aligner = TargetAligner(self.test_reference)
        
        self.assertEqual(len(aligner.targets), 2)
        self.assertIn('target1', aligner.targets)
        self.assertIn('target2', aligner.targets)
        self.assertEqual(len(aligner.targets['target1']), 32)
    
    def test_get_target_regions(self):
        """Test getting target region information"""
        aligner = TargetAligner(self.test_reference)
        regions = aligner.get_target_regions()
        
        self.assertEqual(len(regions), 2)
        self.assertEqual(regions['target1'][1], 32)


class TestVariantCaller(unittest.TestCase):
    """Test variant calling functionality"""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
    
    def test_variant_caller_init(self):
        """Test variant caller initialization"""
        caller = VariantCaller(min_coverage=10, min_variant_frequency=0.2)
        
        self.assertEqual(caller.min_coverage, 10)
        self.assertEqual(caller.min_variant_frequency, 0.2)


class TestCoverageAnalyzer(unittest.TestCase):
    """Test coverage analysis functionality"""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
    
    def test_coverage_analyzer_init(self):
        """Test coverage analyzer initialization"""
        analyzer = CoverageAnalyzer(min_coverage=10)
        
        self.assertEqual(analyzer.min_coverage, 10)


class TestConfig(unittest.TestCase):
    """Test configuration management"""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
    
    def test_default_config(self):
        """Test default configuration"""
        config = Config()
        
        self.assertEqual(config.get('preprocessing', 'min_quality'), 7.0)
        self.assertEqual(config.get('alignment', 'threads'), 4)
    
    def test_save_and_load_config(self):
        """Test saving and loading configuration"""
        config = Config()
        config.set('preprocessing', 'min_quality', 10.0)
        
        config_file = os.path.join(self.temp_dir, 'config.json')
        config.save_config(config_file)
        
        new_config = Config(config_file)
        self.assertEqual(new_config.get('preprocessing', 'min_quality'), 10.0)


if __name__ == '__main__':
    unittest.main()
