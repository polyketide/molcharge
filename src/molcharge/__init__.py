"""molcharge — net charge versus pH, pKa sites and pI of a small molecule from a SMILES (local, offline)."""
__version__ = "0.1.0"

from .charge import compute, detect_sites_dimorphite  # noqa: E402,F401
from .properties import properties as mol_properties  # noqa: E402,F401  (not `properties`: that would shadow the submodule)
