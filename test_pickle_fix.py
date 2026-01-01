import os
import json
import pickle
import numpy as np
import torch
from collections import namedtuple
import sys
import types

from .libconfig import DATA_HOME
from .residue_constants_base import *

# Fix for pickle module reference issues
# Create a fake residue_constants_base module for pickle loading
class PickleModuleFixer:
    @staticmethod
    def find_module(name, path=None):
        if name == 'residue_constants_base':
            # Return a dummy module that will satisfy pickle's import
            fake_module = types.ModuleType(name)
            # Copy all exported items from the real residue_constants_base
            for item in dir():
                if not item.startswith('_'):
                    try:
                        setattr(fake_module, item, eval(item))
                    except (NameError, TypeError):
                        pass
            return fake_module
        return None

# Install the module finder before loading pickle
sys.meta_path.insert(0, PickleModuleFixer())

residue_constants_pkl_fn = DATA_HOME / "residue_constants.pkl"
if residue_constants_pkl_fn.exists():
    use_compiled = True
    with open(residue_constants_pkl_fn, "rb") as fp:
        data_dict = pickle.load(fp)
else:
    use_compiled = False
    data_dict = {}

# Remove the fixer after loading
if PickleModuleFixer in sys.meta_path:
    sys.meta_path.remove(PickleModuleFixer())

if use_compiled:
    ATOM_NAME_ALT_s = data_dict["ATOM_NAME_ALT_s"]
    #
    torsion_s = data_dict["torsion_s"]
else:
    # Fallback to reading from files if pickle not available
    ATOM_NAME_ALT_s = {}
    # ... (file reading code)





