import sys
import json
import pandas as pd
from pathlib import Path
from itertools import product

sys.path.append(str(Path(__file__).parent))
from protein_diversity_pipeline import ProteinDiversityPipeline

def run_parameter_sweep(input_file, output_base_dir="parameter_sweep"):
    """Run pipeline with different parameter combinations"""
    
    # Parameter grid
    n_select_values = [5, 10, 15]
    methods = ['greedy', 'clustering']
    weight_combinations = [
        {'sequence': 0.6, 'physicochemical': 0.2, 'functional': 0.2},
        {'sequence': 0.4, 'physicochemical': 0.3, 'functional': 0.3},
        {'sequence': 0.2, 'physicochemical': 0.4, 'functional': 0.4},
        {'sequence': 0.33, 'physicochemical': 0.33, 'functional': 0.34}
    ]
    
    results_summary = []
    base_path = Path(output_base_dir)
    base_path.mkdir(exist_ok=True)
    
    total_runs = len(n_select_values) * len(methods) * len(weight_combinations)
    run_count = 0
    
    print(f"Running parameter sweep with {total_runs} combinations...")
    
    for n_select, method, weights in product(n_select_values, methods, weight_combinations):
        run_count += 1
        
        # Create unique output directory
        weight_str = f"seq{weights['sequence']:.1f}_phys{weights['physicochemical']:.1f}_func{weights['functional']:.1f}"
        output_dir = base_path / f"run_{run_count:03d}_{method}_n{n_select}_{weight_str}"
        
        print(f"\nRun {run_count}/{total_runs}: {method}, n={n_select}, weights={weight_str}")
        
        try:
            # Initialize pipeline
            pipeline = ProteinDiversityPipeline(str(output_dir))
            pipeline.diversity_optimizer.weights = weights
            
            # Run pipeline
            results = pipeline.run_pipeline(
                input_file=input_file,
                n_select=n_select,
                method=method
            )
            
            # Store summary
            summary = {
                'run_id': run_count,
                'method': method,
                'n_select': n_select,
                'weights': weights,
                'output_dir': str(output_dir),
                'mean_diversity': results['diversity_analysis']['mean_pairwise_diversity'],
                'percentile_rank': results['validation']['percentile_rank'],
                'z_score': results['validation']['z_score'],
                'total_proteins': results['input_info']['total_proteins'],
                'functional_diversity': results['functional_analysis']['unique_ec_numbers'] + 
                                      results['functional_analysis']['unique_go_terms'],
                'success': True
            }
            
        except Exception as e:
            print(f"  Error: {e}")
            summary = {
                'run_id': run_count,
                'method': method,
                'n_select': n_select,
                'weights': weights,
                'output_dir': str(output_dir),
                'error': str(e),
                'success': False
            }
        
        results_summary.append(summary)
    
    # Save summary results
    summary_df = pd.DataFrame(results_summary)
    summary_df.to_csv(base_path / "parameter_sweep_summary.csv", index=False)
    
    with open(base_path / "parameter_sweep_summary.json", "w") as f:
        json.dump(results_summary, f, indent=2)
    
    print(f"\nParameter sweep completed. Summary saved to {base_path}")
    return results_summary

def analyze_parameter_sweep_results(summary_file="parameter_sweep/parameter_sweep_summary.csv"):
    """Analyze results from parameter sweep"""
    
    df = pd.read_csv(summary_file)
    successful_runs = df[df['success'] == True]
    
    if len(successful_runs) == 0:
        print("No successful runs found")
        return
    
    print("Parameter Sweep Analysis")
    print("=" * 30)
    print(f"Total runs: {len(df)}")
    print(f"Successful runs: {len(successful_runs)}")
    print(f"Failed runs: {len(df) - len(successful_runs)}")
    
    # Best performing configurations
    print("\nTop 5 configurations by percentile rank:")
    top_configs = successful_runs.nlargest(5, 'percentile_rank')
    for _, row in top_configs.iterrows():
        print(f"  Run {row['run_id']:3d}: {row['method']:10s} n={row['n_select']:2d} "
              f"rank={row['percentile_rank']:5.1f}% z={row['z_score']:5.2f}")
    
    # Method comparison
    print("\nMethod comparison:")
    method_stats = successful_runs.groupby('method').agg({
        'percentile_rank': ['mean', 'std', 'count'],
        'z_score': ['mean', 'std']
    }).round(2)
    print(method_stats)
    
    # Selection size impact
    print("\nSelection size impact:")
    size_stats = successful_runs.groupby('n_select').agg({
        'percentile_rank': ['mean', 'std'],
        'mean_diversity': ['mean', 'std']
    }).round(3)
    print(size_stats)
    
    return successful_runs
