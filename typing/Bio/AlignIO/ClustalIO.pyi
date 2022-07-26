"""
This type stub file was generated by pyright.
"""

from Bio.AlignIO.Interfaces import AlignmentIterator, SequentialAlignmentWriter

"""Bio.AlignIO support for "clustal" output from CLUSTAL W and other tools.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).
"""
class ClustalWriter(SequentialAlignmentWriter):
    """Clustalw alignment writer."""
    def write_alignment(self, alignment):
        """Use this to write (another) single alignment to an open file."""
        ...
    


class ClustalIterator(AlignmentIterator):
    """Clustalw alignment iterator."""
    _header = ...
    def __next__(self): # -> MultipleSeqAlignment:
        """Parse the next alignment from the handle."""
        ...
    

