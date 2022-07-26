"""
This type stub file was generated by pyright.
"""

from .Interfaces import SequenceIterator

"""Bio.SeqIO support for the "gck" file format.

The GCK binary format is generated by the Gene Construction Kit software
from Textco BioSoftware, Inc.
"""
class GckIterator(SequenceIterator):
    """Parser for GCK files."""
    def __init__(self, source) -> None:
        """Break up a GCK file into SeqRecord objects."""
        ...
    
    def parse(self, handle): # -> Generator[SeqRecord, None, None]:
        """Start parsing the file, and return a SeqRecord generator.

        Note that a GCK file can only contain one sequence, so this
        iterator will always return a single record.
        """
        ...
    

