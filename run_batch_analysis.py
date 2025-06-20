#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent))

from batch_analysis import run_parameter_sweep, analyze_parameter_sweep_results

def main():
    input_file = sys.argv[1] if len(sys.argv) > 1 else "proteins.fasta"
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "batch_results"
    
    print(f"Running batch analysis on {input_file}")
    print("This will test multiple parameter combinations...")
    
    # Run parameter sweep
    results = run_parameter_sweep(input_file, output_dir)
    
    # Analyze results
    print("\n" + "="*50)
    print("BATCH ANALYSIS COMPLETE")
    print("="*50)
    
    analysis = analyze_parameter_sweep_results(f"{output_dir}/parameter_sweep_summary.csv")
    
    print(f"\nDetailed results saved to: {output_dir}/")

if __name__ == "__main__":
    main()
