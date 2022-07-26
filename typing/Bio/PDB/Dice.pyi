"""
This type stub file was generated by pyright.
"""

"""Code for chopping up (dicing) a structure.

This module is used internally by the Bio.PDB.extract() function.
"""
_hydrogen = ...
class ChainSelector:
    """Only accepts residues with right chainid, between start and end.

    Remove hydrogens, waters and ligands. Only use model 0 by default.
    """
    def __init__(self, chain_id, start, end, model_id=...) -> None:
        """Initialize the class."""
        ...
    
    def accept_model(self, model): # -> Literal[1, 0]:
        """Verify if model match the model identifier."""
        ...
    
    def accept_chain(self, chain): # -> Literal[1, 0]:
        """Verify if chain match chain identifier."""
        ...
    
    def accept_residue(self, residue): # -> Literal[0, 1]:
        """Verify if a residue sequence is between the start and end sequence."""
        ...
    
    def accept_atom(self, atom): # -> Literal[0, 1]:
        """Verify if atoms are not Hydrogen."""
        ...
    


def extract(structure, chain_id, start, end, filename): # -> None:
    """Write out selected portion to filename."""
    ...
