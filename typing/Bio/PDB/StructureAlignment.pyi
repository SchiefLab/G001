"""
This type stub file was generated by pyright.
"""

"""Map residues of two structures to each other based on a FASTA alignment."""
class StructureAlignment:
    """Class to align two structures based on an alignment of their sequences."""
    def __init__(self, fasta_align, m1, m2, si=..., sj=...) -> None:
        """Initialize.

        Attributes:
         - fasta_align - Alignment object
         - m1, m2 - two models
         - si, sj - the sequences in the Alignment object that
           correspond to the structures

        """
        ...
    
    def get_maps(self): # -> tuple[dict[Unknown, Unknown], dict[Unknown, Unknown]]:
        """Map residues between the structures.

        Return two dictionaries that map a residue in one structure to
        the equivealent residue in the other structure.
        """
        ...
    
    def get_iterator(self): # -> Generator[Unknown, None, None]:
        """Create an iterator over all residue pairs."""
        ...
    


