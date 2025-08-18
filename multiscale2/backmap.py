import os
import glob
import numpy as np
import MDAnalysis as mda
from .utils import run_command
import yaml
import itertools

def get_latest_pdb(directory):
    """Finds the CALVADOS final structure file with timestamp pattern."""
    # CALVADOS final structure pattern: {sysname}_{YYYYMMDD_HHhMMmSSs}.pdb
    timestamp_pattern = '*_*_*h*m*s.pdb'
    timestamp_files = glob.glob(os.path.join(directory, timestamp_pattern))
    
    if timestamp_files:
        # Find the latest timestamp file
        latest_file = max(timestamp_files, key=os.path.getctime)
        print(f"Found CALVADOS final structure: {os.path.basename(latest_file)}")
        return latest_file
    
    print(f"No CALVADOS final structure found in directory: {directory}")
    return None

def get_protein_dimensions(universe):
    """Calculates the dimensions of the protein atoms."""
    if not universe.atoms:
        return None
    positions = universe.atoms.positions
    min_coords = positions.min(axis=0)
    max_coords = positions.max(axis=0)
    return max_coords - min_coords

def check_protein_fits_in_box(protein_dims, box_dims):
    """Checks if protein dimensions fit within the box dimensions."""
    return np.all(protein_dims <= box_dims)

def write_pdb_with_bfactors(universe, output_pdb_path):
    """
    Writes a PDB file from an MDAnalysis Universe, ensuring the B-factor column is present.
    """
    # Ensure b-factors exist, setting to 0.0 if not.
    if not hasattr(universe.atoms, 'bfactors'):
        universe.atoms.bfactors = np.zeros(universe.atoms.n_atoms)

    with mda.Writer(output_pdb_path, multiframe=False) as W:
        W.write(universe.atoms)

def resize_pdb_box(input_pdb, output_pdb, new_box_dims_A):
    """
    Resizes the box dimensions using MDAnalysis, verifies protein dimensions,
    and writes a new PDB file with updated CRYST1 record.
    """
    # Use MDAnalysis for robust dimension checking
    universe = mda.Universe(input_pdb)
    protein_dims = get_protein_dimensions(universe)
    if protein_dims is None:
        raise ValueError(f"No atoms found in {input_pdb}")

    print(f"Protein dimensions (Angstroms): {protein_dims}")
    print(f"New box dimensions (Angstroms): {new_box_dims_A}")

    if not check_protein_fits_in_box(protein_dims, new_box_dims_A):
        raise ValueError("Protein dimensions exceed the new box dimensions.")

    # Update universe dimensions (in nm, convert from Angstroms)
    new_box_dims_nm = np.array(new_box_dims_A) / 10.0
    universe.dimensions = list(new_box_dims_nm) + [90, 90, 90]  # Add angles
    
    # Write PDB with updated dimensions
    write_pdb_with_bfactors(universe, output_pdb)
    
    print(f"✓ Box resized and protein fits. New PDB saved to: {output_pdb}")
    return output_pdb

def center_condensate(input_pdb, output_pdb):
    """
    Centers the protein assembly in the box and wraps coordinates.
    """
    print(f"--- Centering condensate for {os.path.basename(input_pdb)} ---")
    u = mda.Universe(input_pdb)
    
    protein = u.select_atoms("protein")
    if len(protein) == 0:
        print("⚠️  No protein atoms found, skipping centering.")
        return input_pdb

    box_center = u.dimensions[:3] / 2.0
    protein_cog = protein.center_of_geometry()
    
    translation_vector = box_center - protein_cog
    u.atoms.translate(translation_vector)
    u.atoms.wrap()
    
    write_pdb_with_bfactors(u, output_pdb)
    print(f"✓ Condensate centered. New PDB saved to: {output_pdb}")
    return output_pdb

def _read_fasta_sequences(fasta_path):
    """Minimal FASTA reader: returns dict{name->sequence} (headers sans '>')."""
    sequences = {}
    if not os.path.exists(fasta_path):
        return sequences
    name = None
    seq_chunks = []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    sequences[name] = ''.join(seq_chunks)
                name = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if name is not None:
        sequences[name] = ''.join(seq_chunks)
    return sequences


def _infer_monomer_length(components_data, components_dir):
    """
    Try to infer monomer residue length from FASTA (IDP) or monomer PDB (MDP).
    Returns (monomer_len, source_str).
    """
    # Prefer FASTA when available (IDP)
    ffasta = components_data.get('defaults', {}).get('ffasta')
    system_section = components_data.get('system', {})
    system_name = None
    if isinstance(system_section, dict) and len(system_section) == 1:
        system_name = next(iter(system_section.keys()))

    if ffasta and os.path.exists(ffasta):
        seqs = _read_fasta_sequences(ffasta)
        monomer_len = None
        if system_name and system_name in seqs:
            monomer_len = len(seqs[system_name])
            return monomer_len, f"FASTA:{system_name}"
        if len(seqs) == 1:
            monomer_len = len(next(iter(seqs.values())))
            return monomer_len, "FASTA:<single>"
        if seqs:
            # Fallback to first sequence if multiple
            first_name, first_seq = next(iter(seqs.items()))
            return len(first_seq), f"FASTA:{first_name}"

    # Otherwise try to locate a monomer PDB in an 'input' folder (MDP)
    input_dir = os.path.join(os.path.dirname(components_dir), 'input')
    candidate_paths = []
    if system_name:
        candidate_paths.append(os.path.join(input_dir, f"{system_name}.pdb"))
    # Add any pdb in input dir as fallback
    if os.path.isdir(input_dir):
        for fname in os.listdir(input_dir):
            if fname.lower().endswith('.pdb'):
                candidate_paths.append(os.path.join(input_dir, fname))
    for pdb_path in candidate_paths:
        if os.path.exists(pdb_path):
            u_mono = mda.Universe(pdb_path)
            mono_len = len(u_mono.select_atoms('protein').residues)
            if mono_len > 0:
                return mono_len, f"PDB:{os.path.basename(pdb_path)}"

    return None, None


def add_chain_ids(input_pdb, output_pdb, components_yaml_path):
    """
    Assigns chain IDs to molecules based on a components.yaml file.
    Supports CALVADOS components.yaml structures with either:
    - system: { PROTNAME: { nmol: N } }
    - defaults: { nmol: N }
    - components: [ { nmol: N }, ... ]  (legacy)
    Uses monomer length derived from FASTA (preferred for IDP) or monomer PDB (for MDP).
    """
    print(f"--- Adding chain IDs based on {os.path.basename(components_yaml_path)} ---")
    with open(components_yaml_path, 'r') as f:
        components_data = yaml.safe_load(f)

    nmol = None

    # Preferred: system section with a single protein entry
    system_section = components_data.get('system')
    if isinstance(system_section, dict) and len(system_section) >= 1:
        if len(system_section) == 1:
            first_key = next(iter(system_section.keys()))
            entry = system_section[first_key] or {}
            nmol = entry.get('nmol')
        else:
            raise ValueError("Multiple components in 'system' not supported for automatic chain ID assignment.")

    # Fallback: defaults
    if nmol is None:
        defaults = components_data.get('defaults', {})
        nmol = defaults.get('nmol')

    # Legacy fallback: components list
    if nmol is None:
        comps = components_data.get('components')
        if isinstance(comps, list) and len(comps) > 0:
            nmol = comps[0].get('nmol')

    if not nmol:
        raise ValueError("Could not determine 'nmol' from components.yaml (checked 'system', 'defaults', and 'components').")

    components_dir = os.path.dirname(components_yaml_path)
    monomer_len, src = _infer_monomer_length(components_data, components_dir)
    if not monomer_len:
        raise ValueError("Failed to infer monomer residue length from FASTA or monomer PDB.")
    print(f"  Monomer length: {monomer_len} (source={src})")

    u = mda.Universe(input_pdb)
    protein = u.select_atoms("protein")
    n_residues_total = len(protein.residues)

    expected_total = nmol * monomer_len
    if n_residues_total != expected_total:
        raise ValueError(
            f"Total residues in PDB ({n_residues_total}) != nmol * monomer_len ({expected_total})."
            " Please ensure the PDB matches the specified components and monomer length."
        )

    print(f"  Found {nmol} molecules with {monomer_len} residues each.")

    chain_ids = itertools.cycle('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz')

    # Assign chain IDs per monomer in residue order using slicing
    for i in range(nmol):
        start_residue_idx = i * monomer_len
        end_residue_idx = (i + 1) * monomer_len
        
        # Select residues by index slicing, then get their atoms
        molecule_residues = protein.residues[start_residue_idx:end_residue_idx]
        molecule_atoms = molecule_residues.atoms
        
        cid = next(chain_ids)

        # Assign chain IDs to each atom individually
        for atom in molecule_atoms:
            atom.chainID = cid

    write_pdb_with_bfactors(u, output_pdb)
    print(f"✓ Chain IDs added. New PDB saved to: {output_pdb}")
    return output_pdb


class Backmapper:
    """
    A wrapper for backmapping tools like cg2all.
    """
    def __init__(self, config):
        self.config = config
        self.tool = config.get('tool', 'cg2all')
        self.cg_bead = config.get('cg_bead_name', 'CA')
        self.fix = "--fix" if config.get('fix_sidechains', True) else ""

    def reconstruct(self, input_cg_pdb, output_aa_pdb):
        """
        Performs the reconstruction from CG to AA.
        """
        print("\n" + "="*50)
        print(f"Starting backmapping with {self.tool}")
        print(f"  Input CG PDB: {input_cg_pdb}")
        print(f"  Output AA PDB: {output_aa_pdb}")
        print("="*50)

        if self.tool == 'cg2all':
            command = f"convert_cg2all -p {input_cg_pdb} -o {output_aa_pdb} --cg {self.cg_bead} {self.fix}"
            
            try:
                run_command(command)
                print(f"✓ Backmapping successful. All-atom structure saved to: {output_aa_pdb}")
            except Exception as e:
                print(f"❌ Backmapping with cg2all failed.")
                print(f"  Executed command: {command}")
                print(f"  Error: {e}")
                raise
        else:
            raise NotImplementedError(f"Backmapping tool '{self.tool}' is not supported.")

def read_cryst1_dims(pdb_path):
    """Read CRYST1 a,b,c from a PDB using MDAnalysis. Returns (a,b,c) floats or (None,None,None)."""
    if not os.path.exists(pdb_path):
        return None, None, None
    
    try:
        universe = mda.Universe(pdb_path)
        # Get dimensions from universe (in Angstroms)
        dims = universe.dimensions[:3]
        if dims is not None and all(d > 0 for d in dims):
            return dims[0], dims[1], dims[2]
    except Exception:
        pass
    
    return None, None, None


def write_cryst1_only(input_pdb, output_pdb, dims_abc):
    """Write PDB with specified CRYST1 dimensions using MDAnalysis."""
    a, b, c = dims_abc
    
    # Load universe from input PDB
    universe = mda.Universe(input_pdb)
    
    # Update dimensions (convert from Angstroms to nm for MDAnalysis)
    dims_nm = np.array([a, b, c]) / 10.0
    universe.dimensions = list(dims_nm) + [90, 90, 90]  # Add angles
    
    # Write PDB with updated dimensions
    write_pdb_with_bfactors(universe, output_pdb)


def enforce_xy_keep_z(source_pdb, requested_dims):
    """
    Determine final box dims where X,Y are taken from source PDB dimensions,
    and Z is taken from requested_dims[2]. Returns np.array([ax, ay, cz]).
    """
    req = np.array(requested_dims, dtype=float)
    
    # Load universe and get dimensions
    u = mda.Universe(source_pdb)
    dims = u.dimensions[:3]
    
    if dims is not None and dims[0] > 0 and dims[1] > 0:
        # Use universe dimensions (in Angstroms)
        final = np.array([dims[0], dims[1], req[2]], dtype=float)
        return final
    
    # Fallback: compute bounding box extents
    pos = u.atoms.positions
    span = pos.max(axis=0) - pos.min(axis=0)
    return np.array([span[0], span[1], req[2]], dtype=float)