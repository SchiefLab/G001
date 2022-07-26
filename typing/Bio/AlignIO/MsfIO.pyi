"""
This type stub file was generated by pyright.
"""

from .Interfaces import AlignmentIterator

"""Bio.AlignIO support for GCG MSF format.

The file format was produced by the GCG PileUp and and LocalPileUp tools,
and later tools such as T-COFFEE and MUSCLE support it as an optional
output format.

The original GCG tool would write gaps at ends of each sequence which could
be missing data as tildes (``~``), whereas internal gaps were periods (``.``)
instead. This parser replaces both with minus signs (``-``) for consistency
with the rest of ``Bio.AlignIO``.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).
"""
class MsfIterator(AlignmentIterator):
    """GCG MSF alignment iterator."""
    _header = ...
    def __next__(self): # -> MultipleSeqAlignment:
        """Parse the next alignment from the handle."""
        ...
    

