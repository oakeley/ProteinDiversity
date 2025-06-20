# requirements.txt content:
# biopython>=1.81
# requests>=2.28.0
# numpy>=1.21.0
# scipy>=1.9.0
# pandas>=1.5.0
# scikit-learn>=1.1.0
# matplotlib>=3.6.0
# seaborn>=0.12.0
# tqdm>=4.64.0
# biotite>=0.37.0

import os
import sys
import json
import requests
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union
from dataclasses import dataclass
from collections import Counter, defaultdict
import logging
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Biopython imports
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import substitution_matrices
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Scientific computing
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.stats import entropy
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

# Plotting
import matplotlib.pyplot as plt
import seaborn as sns

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class ProteinData:
    """Container for protein information"""
    id: str
    sequence: str
    uniprot_data: Optional[Dict] = None
    structural_features: Optional[Dict] = None
    functional_features: Optional[Dict] = None
    alphafold_url: Optional[str] = None

class DataCollector:
    """Phase 1: Data Collection & Preprocessing"""
    
    def __init__(self, output_dir: str = "protein_data"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.session = requests.Session()
    
    def read_input_file(self, input_file: str) -> List[str]:
        """Read protein IDs or FASTA sequences from input file"""
        proteins = []
        file_path = Path(input_file)
        
        if not file_path.exists():
            raise FileNotFoundError(f"Input file {input_file} not found")
        
        # Check if it's a FASTA file
        try:
            with open(input_file, 'r') as f:
                sequences = list(SeqIO.parse(f, "fasta"))
            if sequences:
                logger.info(f"Detected FASTA format with {len(sequences)} sequences")
                for i, record in enumerate(sequences):
                    protein_id = record.id if record.id else f"protein_{i+1}"
                    proteins.append(ProteinData(
                        id=protein_id,
                        sequence=str(record.seq)
                    ))
                return proteins
        except:
            pass
        
        # Assume it's a list of UniProt IDs
        with open(input_file, 'r') as f:
            uniprot_ids = [line.strip() for line in f if line.strip()]
        
        logger.info(f"Detected UniProt ID format with {len(uniprot_ids)} IDs")
        return [ProteinData(id=uid, sequence="") for uid in uniprot_ids]
    
    def fetch_uniprot_data(self, protein_id: str) -> Dict:
        """Fetch protein data from UniProt"""
        try:
            url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"
            response = self.session.get(url, timeout=30)
            
            if response.status_code == 200:
                data = response.json()
                return self._parse_uniprot_response(data)
            else:
                logger.warning(f"Failed to fetch data for {protein_id}: {response.status_code}")
                return {}
        except Exception as e:
            logger.error(f"Error fetching data for {protein_id}: {e}")
            return {}
    
    def _parse_uniprot_response(self, data: Dict) -> Dict:
        """Parse UniProt JSON response"""
        parsed = {
            'sequence': data.get('sequence', {}).get('value', ''),
            'length': data.get('sequence', {}).get('length', 0),
            'ec_numbers': [],
            'go_terms': [],
            'keywords': [],
            'organism': data.get('organism', {}).get('scientificName', ''),
            'function': '',
            'catalytic_activity': []
        }
        
        # Extract EC numbers
        if 'proteinDescription' in data:
            for component in data['proteinDescription'].get('recommendedName', {}).get('ecNumbers', []):
                parsed['ec_numbers'].append(component.get('value', ''))
        
        # Extract GO terms
        for go_term in data.get('uniProtKBCrossReferences', []):
            if go_term.get('database') == 'GO':
                parsed['go_terms'].append({
                    'id': go_term['id'],
                    'category': go_term.get('properties', [{}])[0].get('value', ''),
                    'term': go_term.get('properties', [{}])[1].get('value', '') if len(go_term.get('properties', [])) > 1 else ''
                })
        
        # Extract keywords
        for keyword in data.get('keywords', []):
            parsed['keywords'].append(keyword.get('name', ''))
        
        # Extract catalytic activity
        for comment in data.get('comments', []):
            if comment.get('commentType') == 'CATALYTIC ACTIVITY':
                parsed['catalytic_activity'].append(comment.get('reaction', {}).get('name', ''))
            elif comment.get('commentType') == 'FUNCTION':
                parsed['function'] = comment.get('texts', [{}])[0].get('value', '')
        
        return parsed
    
    def get_alphafold_url(self, protein_id: str) -> Optional[str]:
        """Get AlphaFold structure URL if available"""
        # AlphaFold DB follows pattern: AF-{UNIPROT_ID}-F1-model_v4.pdb
        alphafold_url = f"https://alphafold.ebi.ac.uk/files/AF-{protein_id}-F1-model_v4.pdb"
        
        try:
            response = self.session.head(alphafold_url, timeout=10)
            if response.status_code == 200:
                return alphafold_url
        except:
            pass
        return None
    
    def collect_all_data(self, proteins: List[ProteinData]) -> List[ProteinData]:
        """Collect all data for proteins"""
        logger.info("Collecting protein data...")
        
        for protein in tqdm(proteins, desc="Fetching UniProt data"):
            if not protein.sequence:  # Only fetch if we don't have sequence
                uniprot_data = self.fetch_uniprot_data(protein.id)
                protein.uniprot_data = uniprot_data
                protein.sequence = uniprot_data.get('sequence', '')
            
            # Get AlphaFold URL
            protein.alphafold_url = self.get_alphafold_url(protein.id)
        
        # Filter out proteins without sequences
        valid_proteins = [p for p in proteins if p.sequence]
        logger.info(f"Collected data for {len(valid_proteins)} proteins with sequences")
        
        return valid_proteins

class FeatureCalculator:
    """Phase 2: Multi-Dimensional Feature Engineering"""
    
    def __init__(self):
        self.blosum62 = substitution_matrices.load("BLOSUM62")
        self.aligner = Align.PairwiseAligner()
        self.aligner.substitution_matrix = self.blosum62
        self.aligner.open_gap_score = -10
        self.aligner.extend_gap_score = -0.5
    
    def calculate_sequence_features(self, proteins: List[ProteinData]) -> np.ndarray:
        """Calculate sequence-based diversity matrix"""
        logger.info("Calculating sequence features...")
        
        n_proteins = len(proteins)
        diversity_matrix = np.zeros((n_proteins, n_proteins))
        
        # Calculate pairwise BLOSUM-based diversity
        for i in tqdm(range(n_proteins), desc="Sequence alignment"):
            for j in range(i, n_proteins):
                if i == j:
                    diversity_matrix[i, j] = 0
                else:
                    alignments = self.aligner.align(proteins[i].sequence, proteins[j].sequence)
                    score = alignments[0].score if alignments else 0
                    
                    # Convert similarity to diversity (higher score = more similar = less diverse)
                    max_possible = min(len(proteins[i].sequence), len(proteins[j].sequence)) * 11  # Max BLOSUM score
                    diversity = max_possible - score
                    diversity_matrix[i, j] = diversity_matrix[j, i] = max(0, diversity)
        
        return diversity_matrix
    
    def calculate_physicochemical_features(self, proteins: List[ProteinData]) -> np.ndarray:
        """Calculate physicochemical property diversity"""
        logger.info("Calculating physicochemical features...")
        
        properties = []
        for protein in proteins:
            try:
                analyzer = ProteinAnalysis(protein.sequence)
                props = {
                    'molecular_weight': analyzer.molecular_weight(),
                    'isoelectric_point': analyzer.isoelectric_point(),
                    'instability_index': analyzer.instability_index(),
                    'gravy': analyzer.gravy(),  # Grand average of hydropathy
                    'aromaticity': analyzer.aromaticity(),
                    'secondary_structure_fraction': analyzer.secondary_structure_fraction()
                }
                # Flatten secondary structure
                props.update({
                    'helix_fraction': props['secondary_structure_fraction'][0],
                    'turn_fraction': props['secondary_structure_fraction'][1],
                    'sheet_fraction': props['secondary_structure_fraction'][2]
                })
                del props['secondary_structure_fraction']
                properties.append(props)
            except Exception as e:
                logger.warning(f"Error calculating properties for {protein.id}: {e}")
                # Use default values
                properties.append({
                    'molecular_weight': 50000,
                    'isoelectric_point': 7.0,
                    'instability_index': 40.0,
                    'gravy': 0.0,
                    'aromaticity': 0.1,
                    'helix_fraction': 0.3,
                    'turn_fraction': 0.3,
                    'sheet_fraction': 0.4
                })
        
        # Convert to matrix and calculate pairwise distances
        prop_matrix = np.array([[p[key] for key in properties[0].keys()] for p in properties])
        
        # Standardize features
        scaler = StandardScaler()
        prop_matrix_scaled = scaler.fit_transform(prop_matrix)
        
        # Calculate pairwise Euclidean distances
        distances = pdist(prop_matrix_scaled, metric='euclidean')
        diversity_matrix = squareform(distances)
        
        return diversity_matrix
    
    def calculate_functional_features(self, proteins: List[ProteinData]) -> np.ndarray:
        """Calculate functional diversity based on annotations"""
        logger.info("Calculating functional features...")
        
        n_proteins = len(proteins)
        diversity_matrix = np.zeros((n_proteins, n_proteins))
        
        for i in range(n_proteins):
            for j in range(i, n_proteins):
                if i == j:
                    diversity_matrix[i, j] = 0
                else:
                    div_score = self._calculate_functional_diversity(proteins[i], proteins[j])
                    diversity_matrix[i, j] = diversity_matrix[j, i] = div_score
        
        return diversity_matrix
    
    def _calculate_functional_diversity(self, protein1: ProteinData, protein2: ProteinData) -> float:
        """Calculate functional diversity between two proteins"""
        diversity_score = 0.0
        
        # EC number diversity
        if (protein1.uniprot_data and protein2.uniprot_data and 
            protein1.uniprot_data.get('ec_numbers') and protein2.uniprot_data.get('ec_numbers')):
            
            ec1 = protein1.uniprot_data['ec_numbers']
            ec2 = protein2.uniprot_data['ec_numbers']
            ec_diversity = self._calculate_ec_diversity(ec1, ec2)
            diversity_score += ec_diversity * 100  # Weight EC differences highly
        
        # GO term diversity
        if (protein1.uniprot_data and protein2.uniprot_data and 
            protein1.uniprot_data.get('go_terms') and protein2.uniprot_data.get('go_terms')):
            
            go1 = set(term['id'] for term in protein1.uniprot_data['go_terms'])
            go2 = set(term['id'] for term in protein2.uniprot_data['go_terms'])
            
            if go1 and go2:
                jaccard_similarity = len(go1.intersection(go2)) / len(go1.union(go2))
                go_diversity = 1 - jaccard_similarity
                diversity_score += go_diversity * 50  # Moderate weight for GO terms
        
        # Keyword diversity
        if (protein1.uniprot_data and protein2.uniprot_data and 
            protein1.uniprot_data.get('keywords') and protein2.uniprot_data.get('keywords')):
            
            kw1 = set(protein1.uniprot_data['keywords'])
            kw2 = set(protein2.uniprot_data['keywords'])
            
            if kw1 and kw2:
                jaccard_similarity = len(kw1.intersection(kw2)) / len(kw1.union(kw2))
                kw_diversity = 1 - jaccard_similarity
                diversity_score += kw_diversity * 25  # Lower weight for keywords
        
        return diversity_score
    
    def _calculate_ec_diversity(self, ec_list1: List[str], ec_list2: List[str]) -> float:
        """Calculate EC number diversity with hierarchical weighting"""
        if not ec_list1 or not ec_list2:
            return 1.0
        
        max_diversity = 0.0
        weights = [1000, 100, 10, 1]  # Weights for EC hierarchy levels
        
        for ec1 in ec_list1:
            for ec2 in ec_list2:
                ec1_parts = ec1.split('.')
                ec2_parts = ec2.split('.')
                
                diversity = 0
                for i, (part1, part2) in enumerate(zip(ec1_parts, ec2_parts)):
                    if part1 != part2:
                        diversity += weights[min(i, len(weights)-1)]
                
                max_diversity = max(max_diversity, diversity)
        
        return min(max_diversity / 1000, 1.0)  # Normalize to 0-1

class DiversityOptimizer:
    """Phase 3 & 4: Integrated Diversity Scoring and Selection"""
    
    def __init__(self, weights: Dict[str, float] = None):
        self.weights = weights or {'sequence': 0.4, 'physicochemical': 0.3, 'functional': 0.3}
    
    def integrate_diversity_matrices(self, matrices: Dict[str, np.ndarray]) -> np.ndarray:
        """Combine multiple diversity matrices"""
        logger.info("Integrating diversity matrices...")
        
        # Normalize each matrix to 0-1 scale
        normalized_matrices = {}
        for name, matrix in matrices.items():
            if matrix.max() > 0:
                normalized_matrices[name] = matrix / matrix.max()
            else:
                normalized_matrices[name] = matrix
        
        # Weighted combination
        combined_matrix = np.zeros_like(list(normalized_matrices.values())[0])
        for name, matrix in normalized_matrices.items():
            weight = self.weights.get(name, 0.0)
            combined_matrix += weight * matrix
        
        return combined_matrix
    
    def select_diverse_proteins_greedy(self, diversity_matrix: np.ndarray, n_select: int = 10) -> List[int]:
        """Greedy selection for maximum diversity"""
        logger.info(f"Selecting {n_select} most diverse proteins using greedy algorithm...")
        
        n_proteins = diversity_matrix.shape[0]
        if n_select >= n_proteins:
            return list(range(n_proteins))
        
        selected = []
        remaining = list(range(n_proteins))
        
        # Start with protein that has maximum average diversity to all others
        avg_diversity = np.mean(diversity_matrix, axis=1)
        first_protein = np.argmax(avg_diversity)
        selected.append(first_protein)
        remaining.remove(first_protein)
        
        # Iteratively select protein that maximizes minimum distance to selected set
        for _ in tqdm(range(n_select - 1), desc="Greedy selection"):
            best_candidate = None
            best_min_distance = -1
            
            for candidate in remaining:
                min_distance = min(diversity_matrix[candidate][s] for s in selected)
                if min_distance > best_min_distance:
                    best_min_distance = min_distance
                    best_candidate = candidate
            
            if best_candidate is not None:
                selected.append(best_candidate)
                remaining.remove(best_candidate)
        
        return selected
    
    def select_diverse_proteins_clustering(self, diversity_matrix: np.ndarray, n_select: int = 10) -> List[int]:
        """Clustering-based selection for diversity"""
        logger.info(f"Selecting {n_select} most diverse proteins using clustering...")
        
        n_proteins = diversity_matrix.shape[0]
        if n_select >= n_proteins:
            return list(range(n_proteins))
        
        # Convert diversity to distance matrix for clustering
        distance_matrix = diversity_matrix.copy()
        
        # Hierarchical clustering
        condensed_distances = squareform(distance_matrix)
        linkage_matrix = linkage(condensed_distances, method='ward')
        
        # Cut tree to get desired number of clusters
        clusters = fcluster(linkage_matrix, n_select, criterion='maxclust')
        
        # Select representative from each cluster (most central member)
        selected = []
        for cluster_id in range(1, n_select + 1):
            cluster_members = np.where(clusters == cluster_id)[0]
            
            if len(cluster_members) == 1:
                selected.append(cluster_members[0])
            else:
                # Find most representative member (minimum average distance to cluster)
                cluster_distances = distance_matrix[cluster_members][:, cluster_members]
                avg_distances = np.mean(cluster_distances, axis=1)
                representative_idx = np.argmin(avg_distances)
                selected.append(cluster_members[representative_idx])
        
        return selected

class ValidationAnalyzer:
    """Phase 5: Validation & Analysis"""
    
    def __init__(self, output_dir: str = "results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
    
    def validate_selection(self, diversity_matrix: np.ndarray, selected_indices: List[int], 
                          n_bootstrap: int = 100) -> Dict[str, float]:
        """Bootstrap validation of diversity selection"""
        logger.info(f"Validating selection with {n_bootstrap} bootstrap samples...")
        
        n_proteins = diversity_matrix.shape[0]
        n_selected = len(selected_indices)
        
        # Calculate diversity score for current selection
        selected_matrix = diversity_matrix[np.ix_(selected_indices, selected_indices)]
        current_diversity = np.mean(selected_matrix[np.triu_indices(n_selected, k=1)])
        
        # Bootstrap validation
        bootstrap_scores = []
        for _ in tqdm(range(n_bootstrap), desc="Bootstrap validation"):
            # Random selection of same size
            random_indices = np.random.choice(n_proteins, size=n_selected, replace=False)
            random_matrix = diversity_matrix[np.ix_(random_indices, random_indices)]
            random_diversity = np.mean(random_matrix[np.triu_indices(n_selected, k=1)])
            bootstrap_scores.append(random_diversity)
        
        # Calculate statistics
        bootstrap_mean = np.mean(bootstrap_scores)
        bootstrap_std = np.std(bootstrap_scores)
        percentile_rank = (np.sum(np.array(bootstrap_scores) < current_diversity) / n_bootstrap) * 100
        
        return {
            'current_diversity': float(current_diversity),
            'bootstrap_mean': float(bootstrap_mean),
            'bootstrap_std': float(bootstrap_std),
            'percentile_rank': float(percentile_rank),
            'z_score': float((current_diversity - bootstrap_mean) / bootstrap_std if bootstrap_std > 0 else 0)
        }
    
    def analyze_selection(self, proteins: List[ProteinData], selected_indices: List[int], 
                         diversity_matrix: np.ndarray) -> Dict:
        """Comprehensive analysis of selected proteins"""
        logger.info("Analyzing selected proteins...")
        
        selected_proteins = [proteins[i] for i in selected_indices]
        
        analysis = {
            'selected_proteins': [{'id': p.id, 'index': i} for i, p in zip(selected_indices, selected_proteins)],
            'diversity_stats': self._calculate_diversity_stats(diversity_matrix, selected_indices),
            'sequence_stats': self._calculate_sequence_stats(selected_proteins),
            'functional_analysis': self._analyze_functional_diversity(selected_proteins)
        }
        
        return analysis
    
    def _calculate_diversity_stats(self, diversity_matrix: np.ndarray, selected_indices: List[int]) -> Dict:
        """Calculate diversity statistics for selected proteins"""
        selected_matrix = diversity_matrix[np.ix_(selected_indices, selected_indices)]
        
        # Get upper triangle (excluding diagonal)
        upper_triangle = selected_matrix[np.triu_indices(len(selected_indices), k=1)]
        
        return {
            'mean_pairwise_diversity': float(np.mean(upper_triangle)),
            'std_pairwise_diversity': float(np.std(upper_triangle)),
            'min_pairwise_diversity': float(np.min(upper_triangle)),
            'max_pairwise_diversity': float(np.max(upper_triangle)),
            'median_pairwise_diversity': float(np.median(upper_triangle))
        }
    
    def _calculate_sequence_stats(self, proteins: List[ProteinData]) -> Dict:
        """Calculate sequence statistics"""
        lengths = [len(p.sequence) for p in proteins]
        
        return {
            'length_stats': {
                'mean': float(np.mean(lengths)),
                'std': float(np.std(lengths)),
                'min': int(np.min(lengths)),
                'max': int(np.max(lengths)),
                'range': int(np.max(lengths) - np.min(lengths))
            },
            'sequence_count': len(proteins)
        }
    
    def _analyze_functional_diversity(self, proteins: List[ProteinData]) -> Dict:
        """Analyze functional diversity of selected proteins"""
        all_ec_numbers = []
        all_go_terms = []
        all_keywords = []
        
        for protein in proteins:
            if protein.uniprot_data:
                all_ec_numbers.extend(protein.uniprot_data.get('ec_numbers', []))
                all_go_terms.extend([term['id'] for term in protein.uniprot_data.get('go_terms', [])])
                all_keywords.extend(protein.uniprot_data.get('keywords', []))
        
        return {
            'unique_ec_numbers': len(set(all_ec_numbers)),
            'unique_go_terms': len(set(all_go_terms)),
            'unique_keywords': len(set(all_keywords)),
            'total_annotations': len(all_ec_numbers) + len(all_go_terms) + len(all_keywords)
        }
    
    def create_visualizations(self, diversity_matrix: np.ndarray, selected_indices: List[int], 
                            proteins: List[ProteinData]):
        """Create visualization plots"""
        logger.info("Creating visualizations...")
        
        plt.style.use('default')
        
        # 1. Diversity matrix heatmap
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # Full diversity matrix
        im1 = ax1.imshow(diversity_matrix, cmap='viridis', aspect='auto')
        ax1.set_title('Full Diversity Matrix')
        ax1.set_xlabel('Protein Index')
        ax1.set_ylabel('Protein Index')
        plt.colorbar(im1, ax=ax1)
        
        # Selected proteins highlighted
        selected_matrix = diversity_matrix[np.ix_(selected_indices, selected_indices)]
        im2 = ax2.imshow(selected_matrix, cmap='viridis', aspect='auto')
        ax2.set_title('Selected Proteins Diversity Matrix')
        ax2.set_xlabel('Selected Protein Index')
        ax2.set_ylabel('Selected Protein Index')
        ax2.set_xticks(range(len(selected_indices)))
        ax2.set_yticks(range(len(selected_indices)))
        ax2.set_xticklabels([proteins[i].id for i in selected_indices], rotation=45, ha='right')
        ax2.set_yticklabels([proteins[i].id for i in selected_indices])
        plt.colorbar(im2, ax=ax2)
        
        # Diversity distribution
        upper_triangle = diversity_matrix[np.triu_indices(diversity_matrix.shape[0], k=1)]
        selected_upper = selected_matrix[np.triu_indices(len(selected_indices), k=1)]
        
        ax3.hist(upper_triangle, bins=50, alpha=0.7, label='All pairs', density=True)
        ax3.hist(selected_upper, bins=20, alpha=0.7, label='Selected pairs', density=True)
        ax3.set_xlabel('Diversity Score')
        ax3.set_ylabel('Density')
        ax3.set_title('Diversity Score Distribution')
        ax3.legend()
        
        # Sequence length distribution
        all_lengths = [len(p.sequence) for p in proteins]
        selected_lengths = [len(proteins[i].sequence) for i in selected_indices]
        
        ax4.hist(all_lengths, bins=30, alpha=0.7, label='All proteins', density=True)
        ax4.hist(selected_lengths, bins=10, alpha=0.7, label='Selected proteins', density=True)
        ax4.set_xlabel('Sequence Length')
        ax4.set_ylabel('Density')
        ax4.set_title('Sequence Length Distribution')
        ax4.legend()
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'diversity_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 2. PCA visualization
        try:
            pca = PCA(n_components=2)
            diversity_2d = pca.fit_transform(diversity_matrix)
            
            plt.figure(figsize=(10, 8))
            plt.scatter(diversity_2d[:, 0], diversity_2d[:, 1], alpha=0.6, s=50, c='lightblue', label='All proteins')
            plt.scatter(diversity_2d[selected_indices, 0], diversity_2d[selected_indices, 1], 
                       alpha=0.8, s=100, c='red', label='Selected proteins')
            
            # Annotate selected proteins
            for i, idx in enumerate(selected_indices):
                plt.annotate(proteins[idx].id, (diversity_2d[idx, 0], diversity_2d[idx, 1]), 
                           xytext=(5, 5), textcoords='offset points', fontsize=8)
            
            plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
            plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
            plt.title('PCA of Protein Diversity Space')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(self.output_dir / 'pca_analysis.png', dpi=300, bbox_inches='tight')
            plt.close()
        except Exception as e:
            logger.warning(f"Could not create PCA plot: {e}")

class ProteinDiversityPipeline:
    """Master Pipeline Class"""
    
    def __init__(self, output_dir: str = "protein_diversity_results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Initialize components
        self.data_collector = DataCollector(str(self.output_dir / "data"))
        self.feature_calculator = FeatureCalculator()
        self.diversity_optimizer = DiversityOptimizer()
        self.validator = ValidationAnalyzer(str(self.output_dir))
        
        # Setup logging for pipeline
        log_file = self.output_dir / "pipeline.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    def run_pipeline(self, input_file: str, n_select: int = 10, method: str = 'greedy') -> Dict:
        """Run the complete diversity selection pipeline"""
        logger.info(f"Starting protein diversity selection pipeline")
        logger.info(f"Input file: {input_file}")
        logger.info(f"Selecting {n_select} proteins using {method} method")
        
        try:
            # Phase 1: Data Collection
            logger.info("=== Phase 1: Data Collection ===")
            proteins = self.data_collector.read_input_file(input_file)
            proteins = self.data_collector.collect_all_data(proteins)
            
            if len(proteins) < n_select:
                raise ValueError(f"Not enough proteins with sequences ({len(proteins)}) to select {n_select}")
            
            # Save protein data
            self._save_protein_data(proteins)
            
            # Phase 2: Feature Engineering
            logger.info("=== Phase 2: Feature Engineering ===")
            sequence_diversity = self.feature_calculator.calculate_sequence_features(proteins)
            physicochemical_diversity = self.feature_calculator.calculate_physicochemical_features(proteins)
            functional_diversity = self.feature_calculator.calculate_functional_features(proteins)
            
            # Phase 3: Integration
            logger.info("=== Phase 3: Diversity Integration ===")
            diversity_matrices = {
                'sequence': sequence_diversity,
                'physicochemical': physicochemical_diversity,
                'functional': functional_diversity
            }
            
            combined_diversity = self.diversity_optimizer.integrate_diversity_matrices(diversity_matrices)
            
            # Phase 4: Selection
            logger.info("=== Phase 4: Protein Selection ===")
            if method == 'greedy':
                selected_indices = self.diversity_optimizer.select_diverse_proteins_greedy(
                    combined_diversity, n_select)
            elif method == 'clustering':
                selected_indices = self.diversity_optimizer.select_diverse_proteins_clustering(
                    combined_diversity, n_select)
            else:
                raise ValueError(f"Unknown selection method: {method}")
            
            # Phase 5: Validation and Analysis
            logger.info("=== Phase 5: Validation and Analysis ===")
            validation_results = self.validator.validate_selection(combined_diversity, selected_indices)
            analysis_results = self.validator.analyze_selection(proteins, selected_indices, combined_diversity)
            
            # Create visualizations
            self.validator.create_visualizations(combined_diversity, selected_indices, proteins)
            
            # Compile final results
            final_results = {
                'input_info': {
                    'total_proteins': len(proteins),
                    'selected_count': len(selected_indices),
                    'selection_method': method
                },
                'selected_proteins': analysis_results['selected_proteins'],
                'diversity_analysis': analysis_results['diversity_stats'],
                'sequence_analysis': analysis_results['sequence_stats'],
                'functional_analysis': analysis_results['functional_analysis'],
                'validation': validation_results,
                'diversity_matrices': {
                    'sequence_diversity_shape': sequence_diversity.shape,
                    'physicochemical_diversity_shape': physicochemical_diversity.shape,
                    'functional_diversity_shape': functional_diversity.shape,
                    'combined_diversity_shape': combined_diversity.shape
                }
            }
            
            # Save results
            self._save_results(final_results, proteins, selected_indices, combined_diversity)
            
            logger.info("Pipeline completed successfully!")
            return final_results
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            raise
    
    def _save_protein_data(self, proteins: List[ProteinData]):
        """Save protein data to files"""
        # Save sequences as FASTA
        fasta_file = self.output_dir / "input_proteins.fasta"
        with open(fasta_file, 'w') as f:
            for protein in proteins:
                f.write(f">{protein.id}\n{protein.sequence}\n")
        
        # Save metadata as JSON
        metadata = []
        for protein in proteins:
            meta = {
                'id': protein.id,
                'sequence_length': len(protein.sequence),
                'alphafold_url': protein.alphafold_url,
                'uniprot_data': protein.uniprot_data
            }
            metadata.append(meta)
        
        metadata_file = self.output_dir / "protein_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
    
    def _save_results(self, results: Dict, proteins: List[ProteinData], 
                     selected_indices: List[int], diversity_matrix: np.ndarray):
        """Save all results to files"""
        
        # Convert numpy types to native Python types for JSON serialization
        def convert_numpy_types(obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, dict):
                return {key: convert_numpy_types(value) for key, value in obj.items()}
            elif isinstance(obj, list):
                return [convert_numpy_types(item) for item in obj]
            return obj
        
        # Convert results for JSON serialization
        json_safe_results = convert_numpy_types(results)
        
        # Save main results as JSON
        results_file = self.output_dir / "diversity_selection_results.json"
        with open(results_file, 'w') as f:
            json.dump(json_safe_results, f, indent=2)
        
        # Save selected proteins FASTA
        selected_fasta = self.output_dir / "selected_proteins.fasta"
        with open(selected_fasta, 'w') as f:
            for idx in selected_indices:
                protein = proteins[idx]
                f.write(f">{protein.id}\n{protein.sequence}\n")
        
        # Save diversity matrix
        np.save(self.output_dir / "diversity_matrix.npy", diversity_matrix)
        
        # Save selected indices
        np.save(self.output_dir / "selected_indices.npy", np.array(selected_indices))
        
        # Create summary report
        self._create_summary_report(results, proteins, selected_indices)
    
    def _create_summary_report(self, results: Dict, proteins: List[ProteinData], selected_indices: List[int]):
        """Create a human-readable summary report"""
        report_file = self.output_dir / "summary_report.txt"
        
        with open(report_file, 'w') as f:
            f.write("PROTEIN DIVERSITY SELECTION SUMMARY REPORT\n")
            f.write("=" * 50 + "\n\n")
            
            # Input summary
            f.write(f"INPUT SUMMARY:\n")
            f.write(f"Total proteins analyzed: {results['input_info']['total_proteins']}\n")
            f.write(f"Proteins selected: {results['input_info']['selected_count']}\n")
            f.write(f"Selection method: {results['input_info']['selection_method']}\n\n")
            
            # Selected proteins
            f.write("SELECTED PROTEINS:\n")
            for i, protein_info in enumerate(results['selected_proteins'], 1):
                protein = proteins[protein_info['index']]
                f.write(f"{i:2d}. {protein.id} (Length: {len(protein.sequence)} aa)\n")
                if protein.uniprot_data and protein.uniprot_data.get('function'):
                    f.write(f"    Function: {protein.uniprot_data['function'][:100]}...\n")
                if protein.alphafold_url:
                    f.write(f"    AlphaFold: Available\n")
            f.write("\n")
            
            # Diversity statistics
            f.write("DIVERSITY ANALYSIS:\n")
            div_stats = results['diversity_analysis']
            f.write(f"Mean pairwise diversity: {div_stats['mean_pairwise_diversity']:.3f}\n")
            f.write(f"Std pairwise diversity: {div_stats['std_pairwise_diversity']:.3f}\n")
            f.write(f"Min pairwise diversity: {div_stats['min_pairwise_diversity']:.3f}\n")
            f.write(f"Max pairwise diversity: {div_stats['max_pairwise_diversity']:.3f}\n\n")
            
            # Validation results
            f.write("VALIDATION RESULTS:\n")
            val = results['validation']
            f.write(f"Selection diversity score: {val['current_diversity']:.3f}\n")
            f.write(f"Random selection average: {val['bootstrap_mean']:.3f} ± {val['bootstrap_std']:.3f}\n")
            f.write(f"Percentile rank: {val['percentile_rank']:.1f}% (higher is better)\n")
            f.write(f"Z-score: {val['z_score']:.2f}\n\n")
            
            # Interpretation
            f.write("INTERPRETATION:\n")
            if val['percentile_rank'] > 95:
                f.write("EXCELLENT: Selected proteins are significantly more diverse than random selection.\n")
            elif val['percentile_rank'] > 80:
                f.write("GOOD: Selected proteins show above-average diversity.\n")
            elif val['percentile_rank'] > 50:
                f.write("MODERATE: Selected proteins are somewhat more diverse than random.\n")
            else:
                f.write("POOR: Selected proteins are not significantly more diverse than random selection.\n")
            
            # Functional diversity
            func_stats = results['functional_analysis']
            f.write(f"\nFUNCTIONAL DIVERSITY:\n")
            f.write(f"Unique EC numbers: {func_stats['unique_ec_numbers']}\n")
            f.write(f"Unique GO terms: {func_stats['unique_go_terms']}\n")
            f.write(f"Unique keywords: {func_stats['unique_keywords']}\n")

def main():
    """Main function for command-line usage"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Protein Diversity Selection Pipeline')
    parser.add_argument('input_file', help='Input file (UniProt IDs or FASTA sequences)')
    parser.add_argument('-n', '--n_select', type=int, default=10, 
                       help='Number of proteins to select (default: 10)')
    parser.add_argument('-m', '--method', choices=['greedy', 'clustering'], default='greedy',
                       help='Selection method (default: greedy)')
    parser.add_argument('-o', '--output_dir', default='protein_diversity_results',
                       help='Output directory (default: protein_diversity_results)')
    parser.add_argument('--sequence_weight', type=float, default=0.4,
                       help='Weight for sequence diversity (default: 0.4)')
    parser.add_argument('--physicochemical_weight', type=float, default=0.3,
                       help='Weight for physicochemical diversity (default: 0.3)')
    parser.add_argument('--functional_weight', type=float, default=0.3,
                       help='Weight for functional diversity (default: 0.3)')
    
    args = parser.parse_args()
    
    # Set up custom weights
    weights = {
        'sequence': args.sequence_weight,
        'physicochemical': args.physicochemical_weight,
        'functional': args.functional_weight
    }
    
    # Normalize weights
    total_weight = sum(weights.values())
    weights = {k: v/total_weight for k, v in weights.items()}
    
    # Initialize and run pipeline
    pipeline = ProteinDiversityPipeline(args.output_dir)
    pipeline.diversity_optimizer.weights = weights
    
    try:
        results = pipeline.run_pipeline(
            input_file=args.input_file,
            n_select=args.n_select,
            method=args.method
        )
        
        print(f"\nPipeline completed successfully!")
        print(f"Results saved to: {args.output_dir}")
        print(f"Selected {len(results['selected_proteins'])} proteins:")
        for i, protein in enumerate(results['selected_proteins'], 1):
            print(f"  {i}. {protein['id']}")
        
        print(f"\nValidation score: {results['validation']['percentile_rank']:.1f}th percentile")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
