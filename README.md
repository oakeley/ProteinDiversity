# Protein Diversity Selection Pipeline

A comprehensive pipeline for selecting the most diverse proteins from a dataset based on sequence, structural, and functional features.

## Features

- **Multi-modal analysis**: Combines sequence (BLOSUM), physicochemical, and functional diversity
- **Flexible input**: Accepts UniProt IDs or FASTA sequences
- **Multiple algorithms**: Greedy and clustering-based selection methods
- **Validation**: Bootstrap validation with statistical metrics
- **Visualization**: Automated generation of analysis plots
- **Configurable**: Adjustable weights for different diversity components

## Installation

```bash
# Clone the repository
git clone https://github.com/oakeley/ProteinDiversity.git
cd protein-diversity-pipeline

# Install dependencies
pip install -r requirements.txt

# Or install in development mode
pip install -e .
```

## Quick Start

### Using UniProt IDs
```bash
# Create input file with UniProt IDs (one per line)
echo -e "P00698\nP01308\nP69905\nP68871\nP99999" > my_proteins.txt

# Run pipeline
python protein_diversity_pipeline.py my_proteins.txt -n 3 -m greedy -o results/
```

### Using FASTA sequences
```bash
# Run with FASTA file
python protein_diversity_pipeline.py protein.fasta -n 5 -m clustering
```

### Programmatic usage
```python
from protein_diversity_pipeline import ProteinDiversityPipeline

pipeline = ProteinDiversityPipeline("output_directory")
results = pipeline.protein_diversity_pipeline("input_file.txt", n_select=10, method="greedy")
```

## Command Line Options

```
python protein_diversity_pipeline.py input_file [options]

Required:
  input_file              Input file (UniProt IDs or FASTA sequences)

Optional:
  -n, --n_select         Number of proteins to select (default: 10)
  -m, --method           Selection method: greedy|clustering (default: greedy)
  -o, --output_dir       Output directory (default: protein_diversity_results)
  --sequence_weight      Weight for sequence diversity (default: 0.4)
  --physicochemical_weight Weight for physicochemical diversity (default: 0.3)
  --functional_weight    Weight for functional diversity (default: 0.3)
  --verbose             Verbose output
```

## Output Files

The pipeline generates several output files:

- `diversity_selection_results.json`: Complete results in JSON format
- `selected_proteins.fasta`: Selected protein sequences
- `summary_report.txt`: Human-readable summary
- `diversity_analysis.png`: Visualization plots
- `pca_analysis.png`: PCA plot of diversity space
- `diversity_matrix.npy`: Raw diversity matrix
- `protein_metadata.json`: Input protein metadata

## Methodology

### 1. Data Collection
- Fetches protein data from UniProt
- Retrieves AlphaFold structures when available
- Extracts functional annotations (EC numbers, GO terms)

### 2. Feature Engineering
- **Sequence features**: BLOSUM62-based pairwise alignment scores
- **Physicochemical features**: Molecular weight, pI, hydropathy, etc.
- **Functional features**: EC number similarity, GO term overlap

### 3. Diversity Integration
- Combines multiple diversity matrices with configurable weights
- Normalizes scores to 0-1 scale

### 4. Selection Algorithms
- **Greedy**: Iteratively selects protein with maximum minimum distance to selected set
- **Clustering**: Hierarchical clustering followed by representative selection

### 5. Validation
- Bootstrap analysis comparing to random selection
- Statistical significance testing
- Diversity metrics calculation

## Parameter Tuning

# Complete Guide to Batch Analysis

## What is Batch Analysis?

The `batch_analysis.py` script runs your protein diversity selection pipeline with **multiple parameter combinations** to help you find the optimal settings for your specific dataset. Instead of guessing the best parameters, it systematically tests different approaches and tells you which works best.

## How to Run Batch Analysis

### 1. Basic Usage

```bash
# Run with your FASTA file
python run_batch_analysis.py your_proteins.fasta my_batch_results
```

### 2. Interactive Python Usage

```python
# Start Python interpreter
python

# Import and run
from batch_analysis import run_parameter_sweep, analyze_parameter_sweep_results

# Run the parameter sweep
results = run_parameter_sweep('your_proteins.fasta', 'parameter_sweep_results')

# Analyze results
analysis = analyze_parameter_sweep_results('parameter_sweep_results/parameter_sweep_summary.csv')
```

## What Parameters Are Tested?

The batch analysis tests **36 different combinations**:

### Selection Sizes:
- **5 proteins** (small diverse subset)
- **10 proteins** (medium subset) 
- **15 proteins** (larger subset)

### Methods:
- **Greedy**: Iteratively picks most diverse protein
- **Clustering**: Groups similar proteins, picks representatives

### Weight Combinations:
- **Sequence-heavy** (60% sequence, 20% phys, 20% func): For very similar proteins
- **Balanced** (40% sequence, 30% phys, 30% func): General purpose
- **Function-heavy** (20% sequence, 40% phys, 40% func): When you want catalytic diversity
- **Equal weights** (33% each): Neutral approach

## What Results You Get

### 1. Individual Run Results
Each of the 36 combinations creates its own directory with:
- Selected proteins for that parameter set
- Diversity analysis
- Validation scores
- Visualization plots

### 2. Summary Files

#### `parameter_sweep_summary.csv`
Spreadsheet with all results:
```csv
run_id,method,n_select,weights,mean_diversity,percentile_rank,z_score,success
1,greedy,5,{'sequence':0.6...},2.34,87.3,1.45,True
2,clustering,5,{'sequence':0.6...},2.12,79.1,1.21,True
...
```

#### `parameter_sweep_summary.json`
Same data in JSON format for programmatic analysis

### 3. Analysis Output

When you run the analysis, you get:

```
Parameter Sweep Analysis
==============================
Total runs: 36
Successful runs: 35
Failed runs: 1

Top 5 configurations by percentile rank:
  Run  23: clustering  n=10 rank= 94.2% z= 2.31
  Run  11: greedy      n=15 rank= 91.8% z= 2.14
  Run   7: greedy      n=10 rank= 89.3% z= 1.98
  Run  32: clustering  n= 5 rank= 87.6% z= 1.82
  Run  15: greedy      n= 5 rank= 85.1% z= 1.67

Method comparison:
                percentile_rank                 z_score        
                           mean   std count      mean   std
method                                                      
clustering                82.4  8.3    18      1.34  0.41
greedy                    79.1  9.7    17      1.21  0.48

Selection size impact:
          percentile_rank          mean_diversity        
                     mean   std           mean   std
n_select                                            
5                    78.2  7.8           1.89  0.23
10                   84.7  8.1           2.34  0.31
15                   80.3  9.2           2.78  0.42
```

## How to Interpret Results

### Key Metrics Explained:

1. **Percentile Rank** (most important):
   - **>95%**: Excellent diversity selection
   - **80-95%**: Good diversity selection
   - **50-80%**: Moderate improvement over random
   - **<50%**: Poor selection

2. **Z-Score**:
   - **>2.0**: Statistically significant improvement
   - **1.0-2.0**: Moderate improvement
   - **<1.0**: Weak improvement

3. **Mean Diversity**:
   - Higher = more diverse proteins selected
   - Compare across different n_select values

### Practical Interpretation:

From the example above:
- **Best overall**: Clustering with n=10, 94.2% percentile rank
- **Clustering vs Greedy**: Clustering slightly better (82.4% vs 79.1%)
- **Optimal selection size**: 10 proteins gives best percentile rank
- **Method consistency**: Clustering more consistent (std=8.3 vs 9.7)

## Real-World Example

Let's say you have 60 very similar enzyme variants:

```bash
# Run batch analysis
python run_batch_analysis.py enzyme_variants.fasta enzyme_analysis

# Results might show:
# Best: clustering, n=8, function-heavy weights (20%/40%/40%)
# Percentile rank: 96.8% - your selection is better than 96.8% of random selections
# This means the selected 8 enzymes are extremely diverse catalytically
```

## When to Use Batch Analysis

**Use batch analysis when:**
- You have >20 proteins to choose from
- You're unsure about optimal parameters
- You want to maximize diversity for your specific dataset
- You need to justify your selection scientifically

**Skip batch analysis when:**
- You have <10 proteins (not enough for meaningful comparison)
- You already know your optimal parameters
- You're doing a quick exploratory analysis

## Customizing the Parameter Grid

You can modify the parameters tested by editing `batch_analysis.py`:

```python
# In run_parameter_sweep function, modify these:
n_select_values = [3, 5, 8, 12]  # Custom selection sizes
methods = ['greedy']  # Test only greedy method
weight_combinations = [
    {'sequence': 0.1, 'physicochemical': 0.4, 'functional': 0.5},  # Custom weights
    {'sequence': 0.2, 'physicochemical': 0.3, 'functional': 0.5},
]
```

## Output Directory Structure

After running batch analysis:
```
my_batch_results/
├── parameter_sweep_summary.csv          # Main results table
├── parameter_sweep_summary.json         # Same data as JSON
├── run_001_greedy_n5_seq0.6_phys0.2_func0.2/     # Individual runs
│   ├── diversity_selection_results.json
│   ├── selected_proteins.fasta
│   ├── summary_report.txt
│   └── diversity_analysis.png
├── run_002_clustering_n5_seq0.6_phys0.2_func0.2/
│   └── ...
└── ...  (36 total run directories)
```

## Pro Tips

1. **Start small**: Test with a subset first to estimate runtime
2. **Check failed runs**: Look at error messages in summary
3. **Focus on percentile rank**: It's the most reliable metric
4. **Consider your biology**: Higher functional weights for enzymes, higher sequence weights for very similar proteins
5. **Validate manually**: Check if the "best" selection makes biological sense

This batch analysis approach takes the guesswork out of parameter selection and gives you confidence that you've found the most diverse subset possible!

## Testing

Run the test suite:

```bash
python test_pipeline.py
```

## Configuration

Create a `config.yaml` file to set default parameters:

```yaml
pipeline:
  output_dir: "results"
  n_bootstrap: 100

weights:
  sequence: 0.4
  physicochemical: 0.3
  functional: 0.3

selection:
  method: "greedy"
  n_select: 10
```

## Troubleshooting

### Common Issues

1. **UniProt API errors**: Check internet connection and try again
2. **Memory issues**: Reduce number of input proteins or use clustering method
3. **No functional data**: Some proteins may lack annotations in UniProt

### Performance Tips

- For large datasets (>100 proteins), use clustering method
- Adjust weights based on your specific use case
- Use `--verbose` flag for debugging

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{protein_diversity_pipeline,
  title={Multi-Modal Protein Diversity Selection Pipeline},
  author={Edward J. Oakeley},
  year={2025},
  url={https://github.com/oakeley/ProteinDiversity}
}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.
"""

if __name__ == "__main__":
    print("Supporting scripts and configuration files created!")
    print("\nTo run examples:")
    print("python example_usage.py")
    print("\nTo run the main pipeline:")
    print("python protein_diversity_pipeline.py input_file.txt")
    print("\nTo run tests:")
    print("python test_pipeline.py")# setup.py - Installation script
from setuptools import setup, find_packages

setup(
    name="protein-diversity-selector",
    version="1.0.0",
    description="Multi-modal protein diversity selection pipeline for biocatalysis optimization",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.81",
        "requests>=2.28.0",
        "numpy>=1.21.0",
        "scipy>=1.9.0",
        "pandas>=1.5.0",
        "scikit-learn>=1.1.0",
        "matplotlib>=3.6.0",
        "seaborn>=0.12.0",
        "tqdm>=4.64.0",
        "biotite>=0.37.0"
    ],
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "protein-diversity-select=protein_diversity_pipeline:main",
        ],
    },
)
