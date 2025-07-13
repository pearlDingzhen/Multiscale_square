import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
import numpy as np

def generate_plumed_dat(ref_pdb, domains, num_monomers, output_file):
    """
    Generates a plumed.dat file for contact map restraints.
    This function refactors the logic from contactmap_restrain.py.
    """
    ref_universe = mda.Universe(ref_pdb)
    domain_results = []

    for domain_start, domain_end in domains:
        heavy_atoms = ref_universe.select_atoms(f'not name H* and resid {domain_start}-{domain_end}')
        pairs, distances = capped_distance(
            heavy_atoms.positions,
            heavy_atoms.positions,
            max_cutoff=4.5,
            return_distances=True
        )
        atom1 = heavy_atoms[pairs[:, 0]]
        atom2 = heavy_atoms[pairs[:, 1]]
        mask = (atom1.resids != atom2.resids) & (np.abs(atom1.resids - atom2.resids) > 1)
        
        selected_pairs = pairs[mask]
        selected_distances = distances[mask]
        
        global_indices1 = heavy_atoms[selected_pairs[:, 0]].indices
        global_indices2 = heavy_atoms[selected_pairs[:, 1]].indices
        
        # Remove duplicates
        unique_pairs = set()
        keep_indices = []
        for i, (idx1, idx2) in enumerate(zip(global_indices1, global_indices2)):
            pair = tuple(sorted((idx1, idx2)))
            if pair not in unique_pairs:
                unique_pairs.add(pair)
                keep_indices.append(i)
        
        final_indices = np.column_stack((global_indices1[keep_indices], global_indices2[keep_indices]))
        final_distances = selected_distances[keep_indices]
        domain_results.append(np.column_stack((final_indices, final_distances)))

    num_atoms_in_monomer = len(ref_universe.atoms)
    num_domains = len(domains)

    with open(output_file, 'w') as f:
        for monomer_id in range(num_monomers):
            for domain_id, results in enumerate(domain_results):
                if len(results) == 0: continue
                weight = 1.0 / len(results)
                q_label = monomer_id * num_domains + domain_id
                f.write(f'q{q_label}: CONTACTMAP ...\n')
                
                group_lines = []
                for i, result in enumerate(results):
                    atom1 = int(result[0]) + 1 + monomer_id * num_atoms_in_monomer
                    atom2 = int(result[1]) + 1 + monomer_id * num_atoms_in_monomer
                    ref_value = result[2] / 10.0  # Angstrom to nm
                    
                    line = (f"ATOMS{i+1}={atom1},{atom2} "
                            f"SWITCH{i+1}={{Q R_0=0.01 BETA=20 LAMBDA=1.5 REF={ref_value:.6f}}} "
                            f"WEIGHT{i+1}={weight:.6f}")
                    group_lines.append(line)
                
                f.write(" ".join(group_lines) + "\n")
                f.write('... CONTACTMAP\n')
                f.write(f'PRINT ARG=q{q_label} FILE=COLVAR_{q_label}\n\n')

        f.write("\n# Restraints\n")
        for monomer_id in range(num_monomers * num_domains):
            f.write(f'RESTRAINT ARG=q{monomer_id} AT=1.0 KAPPA=10000 SLOPE=0.\n')