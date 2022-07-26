"""
This type stub file was generated by pyright.
"""

from .Interfaces import SequenceIterator

"""Bio.SeqIO support for the "ig" (IntelliGenetics or MASE) file format.

This module is for reading and writing IntelliGenetics format files as
SeqRecord objects.  This file format appears to be the same as the MASE
multiple sequence alignment format.

You are expected to use this module via the Bio.SeqIO functions.
"""
class IgIterator(SequenceIterator):
    """Parser for IntelliGenetics files."""
    def __init__(self, source) -> None:
        """Iterate over IntelliGenetics records (as SeqRecord objects).

        source - file-like object opened in text mode, or a path to a file

        The optional free format file header lines (which start with two
        semi-colons) are ignored.

        The free format commentary lines at the start of each record (which
        start with a semi-colon) are recorded as a single string with embedded
        new line characters in the SeqRecord's annotations dictionary under the
        key 'comment'.

        Examples
        --------
        >>> with open("IntelliGenetics/TAT_mase_nuc.txt") as handle:
        ...     for record in IgIterator(handle):
        ...         print("%s length %i" % (record.id, len(record)))
        ...
        A_U455 length 303
        B_HXB2R length 306
        C_UG268A length 267
        D_ELI length 309
        F_BZ163A length 309
        O_ANT70 length 342
        O_MVP5180 length 348
        CPZGAB length 309
        CPZANT length 309
        A_ROD length 390
        B_EHOA length 420
        D_MM251 length 390
        STM_STM length 387
        VER_AGM3 length 354
        GRI_AGM677 length 264
        SAB_SAB1C length 219
        SYK_SYK length 330

        """
        ...
    
    def parse(self, handle): # -> Generator[SeqRecord, None, None]:
        """Start parsing the file, and return a SeqRecord generator."""
        ...
    
    def iterate(self, handle): # -> Generator[SeqRecord, None, None]:
        """Iterate over the records in the IntelliGenetics file."""
        ...
    


if __name__ == "__main__":
    ...
