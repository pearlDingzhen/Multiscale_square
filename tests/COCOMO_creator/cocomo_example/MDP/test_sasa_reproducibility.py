#!/usr/bin/env python3
"""
Test SASA calculation reproducibility.
Compare two consecutive SASA calculations to check if they are deterministic.
"""
import numpy as np
from mdsim import PDBReader

def test_sasa_reproducibility(pdb_file, n_runs=5, n_sphere_points=1920):
    """
    Test if SASA calculation is reproducible.
    
    Args:
        pdb_file: Path to PDB file
        n_runs: Number of times to run the calculation
        n_sphere_points: Number of sphere points for SASA calculation
    
    Returns:
        dict: Statistics about reproducibility
    """
    print("=" * 70)
    print("SASA Calculation Reproducibility Test")
    print("=" * 70)
    print(f"PDB file: {pdb_file}")
    print(f"Number of runs: {n_runs}")
    print(f"Sphere points: {n_sphere_points}")
    print()
    
    # Run SASA calculations multiple times
    all_results = []
    
    for run in range(n_runs):
        print(f"Run {run + 1}/{n_runs}...", end=" ", flush=True)
        reader = PDBReader(pdb_file)
        model = reader[0]
        sasa = model.sasa_by_residue(n_sphere_points=n_sphere_points)
        all_results.append(np.array(sasa))
        print(f"Done. Got {len(sasa)} residues.")
    
    print()
    
    # Compare all runs
    print("=" * 70)
    print("Comparison Results")
    print("=" * 70)
    
    # Convert to numpy array for easier comparison
    sasa_array = np.array(all_results)  # Shape: (n_runs, n_residues)
    n_residues = sasa_array.shape[1]
    
    # Calculate statistics
    mean_sasa = np.mean(sasa_array, axis=0)
    std_sasa = np.std(sasa_array, axis=0)
    max_sasa = np.max(sasa_array, axis=0)
    min_sasa = np.min(sasa_array, axis=0)
    
    # Find differences
    max_diff = max_sasa - min_sasa
    max_diff_idx = np.argmax(max_diff)
    max_diff_value = max_diff[max_diff_idx]
    
    # Compare first two runs in detail
    diff_01 = np.abs(sasa_array[0] - sasa_array[1])
    max_diff_01 = np.max(diff_01)
    max_diff_01_idx = np.argmax(diff_01)
    mean_diff_01 = np.mean(diff_01)
    
    print(f"\n1. Overall Statistics:")
    print(f"   Number of residues: {n_residues}")
    print(f"   Mean SASA across all runs: {np.mean(mean_sasa):.6f} nm²")
    print(f"   Mean std across runs: {np.mean(std_sasa):.6f} nm²")
    print(f"   Max std across runs: {np.max(std_sasa):.6f} nm² (residue {np.argmax(std_sasa) + 1})")
    
    print(f"\n2. Run-to-Run Differences:")
    print(f"   Max difference (any two runs): {max_diff_value:.10f} nm²")
    print(f"   Max difference residue: {max_diff_idx + 1}")
    print(f"   Mean difference (run 1 vs run 2): {mean_diff_01:.10f} nm²")
    print(f"   Max difference (run 1 vs run 2): {max_diff_01:.10f} nm² (residue {max_diff_01_idx + 1})")
    
    # Check if all runs are identical
    all_identical = np.allclose(sasa_array[0], sasa_array[1:], atol=1e-10, rtol=1e-10)
    
    print(f"\n3. Reproducibility Check:")
    if all_identical:
        print(f"   ✓ All runs are IDENTICAL (within 1e-10 tolerance)")
        print(f"   → SASA calculation is DETERMINISTIC")
    else:
        print(f"   ✗ Runs show DIFFERENCES")
        print(f"   → SASA calculation may have STOCHASTIC components")
        
        # Show residues with differences
        n_different = np.sum(diff_01 > 1e-10)
        print(f"\n   Residues with differences > 1e-10: {n_different}/{n_residues}")
        
        if n_different > 0:
            print(f"\n   Top 10 residues with largest differences:")
            top_indices = np.argsort(diff_01)[-10:][::-1]
            for idx in top_indices:
                if diff_01[idx] > 1e-10:
                    print(f"     Residue {idx + 1}: diff = {diff_01[idx]:.10f} nm²")
                    print(f"       Run 1: {sasa_array[0, idx]:.10f} nm²")
                    print(f"       Run 2: {sasa_array[1, idx]:.10f} nm²")
    
    # Compare with different sphere point counts
    print(f"\n4. Testing with different sphere point counts:")
    for n_points in [960, 1920, 3840]:
        print(f"   Testing with {n_points} points...", end=" ", flush=True)
        reader = PDBReader(pdb_file)
        model = reader[0]
        sasa_test = model.sasa_by_residue(n_sphere_points=n_points)
        print(f"Done. Mean SASA: {np.mean(sasa_test):.6f} nm²")
    
    print("\n" + "=" * 70)
    
    return {
        'all_results': all_results,
        'mean_sasa': mean_sasa,
        'std_sasa': std_sasa,
        'max_diff': max_diff_value,
        'is_deterministic': all_identical,
        'n_different': n_different if not all_identical else 0
    }


def main():
    # Use the same PDB file as the test
    pdb_file = '/mnt/hdd1/home/tianxj/code/ms2main/Multiscale_square/tests/COCOMO_wrapper/cocomo_example/MDP/system_with_segid.pdb'
    
    # Test reproducibility
    results = test_sasa_reproducibility(pdb_file, n_runs=5, n_sphere_points=1920)
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    if results['is_deterministic']:
        print("✓ SASA calculation is DETERMINISTIC - same input produces same output")
    else:
        print("✗ SASA calculation shows VARIABILITY - may be stochastic or")
        print("  have numerical precision issues")
        print(f"  Max difference observed: {results['max_diff']:.10f} nm²")
        print(f"  Residues with differences: {results['n_different']}")
    print("=" * 70)


if __name__ == '__main__':
    main()

