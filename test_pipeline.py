#!/usr/bin/env python3
"""
Test suite for the protein diversity pipeline
"""

import unittest
import tempfile
import shutil
from pathlib import Path
import numpy as np
import sys

sys.path.append(str(Path(__file__).parent))
from protein_diversity_pipeline import (
    ProteinData, DataCollector, FeatureCalculator, 
    DiversityOptimizer, ValidationAnalyzer, ProteinDiversityPipeline
)

class TestProteinDiversityPipeline(unittest.TestCase):
    
    def setUp(self):
        """Set up test fixtures"""
        self.temp_dir = tempfile.mkdtemp()
        self.test_sequences = [
            ProteinData("test1", "MKTVRQERLKSIVRILPGKWKRK"),
            ProteinData("test2", "MVLSEGEWQLVLHVWAKVEADVA"),
            ProteinData("test3", "MADNRDPASDQMKHWKEQRAAQK"),
            ProteinData("test4", "IVGGYTCGANTVPYQVSLNSGYHFC"),
            ProteinData("test5", "PTIGPSENSAPGQVQTSQTPEP")
        ]
    
    def tearDown(self):
        """Clean up test fixtures"""
        shutil.rmtree(self.temp_dir)
    
    def test_data_collector_fasta_reading(self):
        """Test reading FASTA files"""
        # Create test FASTA file
        fasta_content = ">test1\nMKTVRQERLKSIVRILPGKWKRK\n>test2\nMVLSEGEWQLVLHVWAKVEADVA\n"
        fasta_file = Path(self.temp_dir) / "test.fasta"
        with open(fasta_file, 'w') as f:
            f.write(fasta_content)
        
        collector = DataCollector(self.temp_dir)
        proteins = collector.read_input_file(str(fasta_file))
        
        self.assertEqual(len(proteins), 2)
        self.assertEqual(proteins[0].id, "test1")
        self.assertEqual(proteins[0].sequence, "MKTVRQERLKSIVRILPGKWKRK")
    
    def test_data_collector_uniprot_ids(self):
        """Test reading UniProt ID files"""
        # Create test ID file
        id_content = "P00698\nP01308\n"
        id_file = Path(self.temp_dir) / "test_ids.txt"
        with open(id_file, 'w') as f:
            f.write(id_content)
        
        collector = DataCollector(self.temp_dir)
        proteins = collector.read_input_file(str(id_file))
        
        self.assertEqual(len(proteins), 2)
        self.assertEqual(proteins[0].id, "P00698")
        self.assertEqual(proteins[1].id, "P01308")
    
    def test_feature_calculator_sequence_features(self):
        """Test sequence feature calculation"""
        calculator = FeatureCalculator()
        diversity_matrix = calculator.calculate_sequence_features(self.test_sequences)
        
        # Check matrix properties
        self.assertEqual(diversity_matrix.shape, (5, 5))
        self.assertTrue(np.allclose(diversity_matrix, diversity_matrix.T))  # Symmetric
        self.assertTrue(np.all(np.diag(diversity_matrix) == 0))  # Zero diagonal
        self.assertTrue(np.all(diversity_matrix >= 0))  # Non-negative
    
    def test_feature_calculator_physicochemical_features(self):
        """Test physicochemical feature calculation"""
        calculator = FeatureCalculator()
        diversity_matrix = calculator.calculate_physicochemical_features(self.test_sequences)
        
        # Check matrix properties
        self.assertEqual(diversity_matrix.shape, (5, 5))
        self.assertTrue(np.allclose(diversity_matrix, diversity_matrix.T))  # Symmetric
        self.assertTrue(np.all(diversity_matrix >= 0))  # Non-negative
    
    def test_diversity_optimizer_greedy_selection(self):
        """Test greedy selection algorithm"""
        # Create test diversity matrix
        np.random.seed(42)
        diversity_matrix = np.random.rand(10, 10)
        diversity_matrix = (diversity_matrix + diversity_matrix.T) / 2  # Make symmetric
        np.fill_diagonal(diversity_matrix, 0)  # Zero diagonal
        
        optimizer = DiversityOptimizer()
        selected = optimizer.select_diverse_proteins_greedy(diversity_matrix, n_select=5)
        
        self.assertEqual(len(selected), 5)
        self.assertEqual(len(set(selected)), 5)  # All unique
        self.assertTrue(all(0 <= idx < 10 for idx in selected))  # Valid indices
    
    def test_diversity_optimizer_clustering_selection(self):
        """Test clustering-based selection"""
        # Create test diversity matrix
        np.random.seed(42)
        diversity_matrix = np.random.rand(10, 10)
        diversity_matrix = (diversity_matrix + diversity_matrix.T) / 2  # Make symmetric
        np.fill_diagonal(diversity_matrix, 0)  # Zero diagonal
        
        optimizer = DiversityOptimizer()
        selected = optimizer.select_diverse_proteins_clustering(diversity_matrix, n_select=5)
        
        self.assertEqual(len(selected), 5)
        self.assertEqual(len(set(selected)), 5)  # All unique
        self.assertTrue(all(0 <= idx < 10 for idx in selected))  # Valid indices
    
    def test_diversity_optimizer_matrix_integration(self):
        """Test diversity matrix integration"""
        # Create test matrices
        np.random.seed(42)
        matrices = {
            'sequence': np.random.rand(5, 5),
            'physicochemical': np.random.rand(5, 5),
            'functional': np.random.rand(5, 5)
        }
        
        optimizer = DiversityOptimizer()
        combined = optimizer.integrate_diversity_matrices(matrices)
        
        self.assertEqual(combined.shape, (5, 5))
        self.assertTrue(np.all(combined >= 0))  # Non-negative
        self.assertTrue(np.all(combined <= 1))  # Normalized
    
    def test_validation_analyzer(self):
        """Test validation analysis"""
        # Create test data
        np.random.seed(42)
        diversity_matrix = np.random.rand(20, 20)
        diversity_matrix = (diversity_matrix + diversity_matrix.T) / 2
        np.fill_diagonal(diversity_matrix, 0)
        
        selected_indices = [0, 5, 10, 15]
        
        validator = ValidationAnalyzer(self.temp_dir)
        validation_results = validator.validate_selection(
            diversity_matrix, selected_indices, n_bootstrap=10
        )
        
        # Check required keys
        required_keys = ['current_diversity', 'bootstrap_mean', 'bootstrap_std', 
                        'percentile_rank', 'z_score']
        for key in required_keys:
            self.assertIn(key, validation_results)
        
        # Check value ranges
        self.assertGreaterEqual(validation_results['percentile_rank'], 0)
        self.assertLessEqual(validation_results['percentile_rank'], 100)
    
    def test_end_to_end_pipeline_fasta(self):
        """Test complete pipeline with FASTA input"""
        # Create test FASTA file
        fasta_content = """
>prot1
MKTVRQERLKSIVRILPGKWKRK
>prot2
MVLSEGEWQLVLHVWAKVEADVA
>prot3
MADNRDPASDQMKHWKEQRAAQK
>prot4
IVGGYTCGANTVPYQVSLNSGYHFC
>prot5
PTIGPSENSAPGQVQTSQTPEP
>prot6
MALWMRLLPLLALLALWGPDPAA
>prot7
MKALIVLGLVLLSVTVQGKVFER
>prot8
MVLSPADKTNVKAAWGKVGAHA
        """.strip()
        
        fasta_file = Path(self.temp_dir) / "test.fasta"
        with open(fasta_file, 'w') as f:
            f.write(fasta_content)
        
        # Run pipeline
        pipeline = ProteinDiversityPipeline(self.temp_dir)
        results = pipeline.run_pipeline(str(fasta_file), n_select=4, method='greedy')
        
        # Check results structure
        self.assertIn('selected_proteins', results)
        self.assertIn('diversity_analysis', results)
        self.assertIn('validation', results)
        self.assertEqual(len(results['selected_proteins']), 4)
        
        # Check output files
        expected_files = [
            'diversity_selection_results.json',
            'selected_proteins.fasta',
            'summary_report.txt',
            'diversity_matrix.npy'
        ]
        
        for filename in expected_files:
            file_path = Path(self.temp_dir) / filename
            self.assertTrue(file_path.exists(), f"Missing output file: {filename}")

def run_tests():
    """Run all tests"""
    unittest.main(verbosity=2)

if __name__ == "__main__":
    run_tests()
