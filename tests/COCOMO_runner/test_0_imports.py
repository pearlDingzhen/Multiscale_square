#!/usr/bin/env python3
"""
Test 0: Basic Import and Functionality Tests for COCOMO2 Creator

This test file validates:
1. All imports work correctly
2. Basic class instantiation
3. Topology information building
4. SASA value loading
"""

import os
import sys
import numpy as np
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))


def test_imports():
    """Test that all required modules can be imported."""
    print("=" * 60)
    print("TEST 0.1: Import Tests")
    print("=" * 60)
    
    try:
        from multiscale2.src import CGSimulationConfig, CGSimulator
        print("  [PASS] CGSimulationConfig and CGSimulator imported")
    except ImportError as e:
        print(f"  [FAIL] Failed to import CGSimulationConfig/CGSimulator: {e}")
        return False
    
    try:
        from multiscale2.src.cocomo2_creator import COCOMO
        print("  [PASS] COCOMO class imported")
    except ImportError as e:
        print(f"  [FAIL] Failed to import COCOMO: {e}")
        return False
    
    try:
        from openmm import (
            Platform, XmlSerializer, VerletIntegrator,
        )
        from openmm.app import DCDReporter, StateDataReporter
        print("  [PASS] OpenMM modules imported")
    except ImportError as e:
        print(f"  [FAIL] Failed to import OpenMM modules: {e}")
        return False
    
    try:
        from openmm.app import PDBFile
        print("  [PASS] OpenMM.app.PDBFile imported")
    except ImportError as e:
        print(f"  [FAIL] Failed to import PDBFile: {e}")
        return False
    
    try:
        import parmed
        print("  [PASS] ParmEd imported")
    except ImportError as e:
        print(f"  [FAIL] Failed to import ParmEd: {e}")
        return False
    
    print()
    return True


def test_config_loading():
    """Test loading configuration from YAML file."""
    print("=" * 60)
    print("TEST 0.2: Configuration Loading")
    print("=" * 60)

    from multiscale2.src import CGSimulationConfig
    
    test_dir = Path(__file__).parent
    config_file = test_dir / "config_idp.yaml"
    
    if not config_file.exists():
        print(f"  [FAIL] Config file not found: {config_file}")
        return False
    
    try:
        config = CGSimulationConfig.from_yaml(str(config_file))
        print(f"  [PASS] Configuration loaded from {config_file}")
        print(f"       System name: {config.system_name}")
        print(f"       Box size: {config.box} nm")
        print(f"       Temperature: {config.temperature} K")
        print(f"       Components: {len(config.components)}")
        
        for comp in config.components:
            print(f"         - {comp.name}: {comp.type.value}, nmol={comp.nmol}")
        
        print()
        return True
    except Exception as e:
        print(f"  [FAIL] Failed to load configuration: {e}")
        return False


def test_simulator_creation():
    """Test creating CGSimulator from config."""
    print("=" * 60)
    print("TEST 0.3: Simulator Creation")
    print("=" * 60)

    from multiscale2.src import CGSimulationConfig, CGSimulator
    
    test_dir = Path(__file__).parent
    config_file = test_dir / "config_idp.yaml"
    output_dir = test_dir / "output_test0_simulator"
    
    # Create output directory
    output_dir.mkdir(exist_ok=True)
    
    try:
        config = CGSimulationConfig.from_yaml(str(config_file))
        sim = CGSimulator(config)
        print(f"  [PASS] CGSimulator created")
        print(f"       Output directory: {output_dir}")
        return True
    except Exception as e:
        print(f"  [FAIL] Failed to create simulator: {e}")
        return False


def test_topology_info():
    """Test building topology information."""
    print("=" * 60)
    print("TEST 0.4: Topology Information Building")
    print("=" * 60)

    from multiscale2.src import CGSimulationConfig, CGSimulator
    
    test_dir = Path(__file__).parent
    config_file = test_dir / "config_idp.yaml"
    output_dir = test_dir / "output_test0_topology"
    output_dir.mkdir(exist_ok=True)
    
    try:
        config = CGSimulationConfig.from_yaml(str(config_file))
        sim = CGSimulator(config)
        sim.setup(str(output_dir), overwrite=True)
        
        # Test global sequence
        global_seq = sim.get_global_sequence()
        print(f"  [PASS] Global sequence generated")
        print(f"       Length: {len(global_seq)}")
        print(f"       First 50 chars: {global_seq[:50]}...")
        
        # Test chain IDs
        chain_ids = sim.get_chain_ids()
        print(f"  [PASS] Chain IDs generated")
        print(f"       Length: {len(chain_ids)}")
        print(f"       Chain range: 1-{max(chain_ids)}")
        
        # Test folded domains
        folded = sim.get_folded_domains()
        print(f"  [PASS] Folded domains generated")
        print(f"       Length: {len(folded)}")
        print(f"       Folded regions: {sum(folded)}")
        
        # Test component names
        comp_names = sim._build_component_names()
        print(f"  [PASS] Component names generated")
        print(f"       Length: {len(comp_names)}")
        print(f"       Unique components: {set(comp_names)}")
        
        print()
        return True
    except Exception as e:
        print(f"  [FAIL] Failed to build topology info: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_cocomo_creation():
    """Test creating COCOMO system."""
    print("=" * 60)
    print("TEST 0.5: COCOMO System Creation")
    print("=" * 60)

    from multiscale2.src import CGSimulationConfig, CGSimulator
    from multiscale2.src.cocomo2_creator import COCOMO
    
    test_dir = Path(__file__).parent
    config_file = test_dir / "config_idp.yaml"
    output_dir = test_dir / "output_test0_cocomo"
    output_dir.mkdir(exist_ok=True)
    
    try:
        from multiscale2.src.cocomo2_creator import COCOMO

        config = CGSimulationConfig.from_yaml(str(config_file))
        sim = CGSimulator(config)
        sim.setup(str(output_dir), overwrite=True)
        
        # Build topology info
        topology_info = {
            'global_sequence': sim.get_global_sequence(),
            'chain_ids': sim.get_chain_ids(),
            'folded_domains': sim.get_folded_domains(),
            'component_names': sim._build_component_names(),
            'local_residue_indices': list(range(1, len(sim.get_global_sequence()) + 1)),
        }
        
        # Create dummy positions (all zeros for testing)
        n_residues = len(topology_info['global_sequence'])
        positions = np.zeros((n_residues, 3))
        
        # Create COCOMO instance
        box_size = config.box
        cocomo = COCOMO(
            box_size=box_size,
            topology_info=topology_info,
            positions=positions,
            surf=0.7,
            resources='CPU'  # Use CPU for testing
        )
        
        print(f"  [PASS] COCOMO instance created")
        print(f"       Box size: {box_size} nm")
        print(f"       Residues: {n_residues}")
        print(f"       Chains: {max(topology_info['chain_ids'])}")
        
        print()
        return True
    except Exception as e:
        print(f"  [FAIL] Failed to create COCOMO system: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_sasa_loading():
    """Test SASA value loading."""
    print("=" * 60)
    print("TEST 0.6: SASA Value Loading")
    print("=" * 60)

    from multiscale2.src import CGSimulationConfig, CGSimulator
    
    test_dir = Path(__file__).parent
    config_file = test_dir / "config_idp.yaml"
    output_dir = test_dir / "output_test0_sasa"
    output_dir.mkdir(exist_ok=True)
    
    try:
        config = CGSimulationConfig.from_yaml(str(config_file))
        sim = CGSimulator(config)
        sim.setup(str(output_dir), overwrite=True)
        
        # Test when SASA file doesn't exist (should return None)
        sasa = sim._get_sasa_values()
        print(f"  [PASS] SASA loading function works")
        print(f"       Return value: {type(sasa)}")
        print(f"       (None expected when file doesn't exist)")
        
        print()
        return True
    except Exception as e:
        print(f"  [FAIL] Failed to test SASA loading: {e}")
        return False


def main():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("COCOMO2 Creator - Test 0: Basic Import and Functionality")
    print("=" * 60 + "\n")
    
    tests = [
        ("Imports", test_imports),
        ("Configuration Loading", test_config_loading),
        ("Simulator Creation", test_simulator_creation),
        ("Topology Info", test_topology_info),
        ("COCOMO Creation", test_cocomo_creation),
        ("SASA Loading", test_sasa_loading),
    ]
    
    results = []
    for name, test_func in tests:
        try:
            passed = test_func()
            results.append((name, passed))
        except Exception as e:
            print(f"  [ERROR] Unexpected error in {name}: {e}")
            results.append((name, False))
    
    # Summary
    print("=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    
    passed = sum(1 for _, r in results if r)
    total = len(results)
    
    for name, result in results:
        status = "PASS" if result else "FAIL"
        print(f"  [{status}] {name}")
    
    print()
    print(f"Total: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n[SUCCESS] All basic tests passed!")
        return 0
    else:
        print("\n[FAILURE] Some tests failed.")
        return 1


if __name__ == '__main__':
    exit(main())

