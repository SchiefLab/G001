"""
This type stub file was generated by pyright.
"""

from .Interfaces import SequenceIterator

"""Bio.SeqIO support for the SnapGene file format.

The SnapGene binary format is the native format used by the SnapGene program
from GSL Biotech LLC.
"""
_packet_handlers = ...
class SnapGeneIterator(SequenceIterator):
    """Parser for SnapGene files."""
    def __init__(self, source) -> None:
        """Parse a SnapGene file and return a SeqRecord object.

        Argument source is a file-like object or a path to a file.

        Note that a SnapGene file can only contain one sequence, so this
        iterator will always return a single record.
        """
        ...
    
    def parse(self, handle): # -> Generator[SeqRecord, None, None]:
        """Start parsing the file, and return a SeqRecord generator."""
        ...
    
    def iterate(self, handle): # -> Generator[SeqRecord, None, None]:
        """Iterate over the records in the SnapGene file."""
        ...
    


