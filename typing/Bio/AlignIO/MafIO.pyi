"""
This type stub file was generated by pyright.
"""

from .Interfaces import SequentialAlignmentWriter

"""Bio.AlignIO support for the "maf" multiple alignment format.

The Multiple Alignment Format, described by UCSC, stores a series of
multiple alignments in a single file. It is suitable for whole-genome
to whole-genome alignments, metadata such as source chromosome, start
position, size, and strand can be stored.

See http://genome.ucsc.edu/FAQ/FAQformat.html#format5

You are expected to use this module via the Bio.AlignIO functions(or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).

Coordinates in the MAF format are defined in terms of zero-based start
positions (like Python) and aligning region sizes.

A minimal aligned region of length one and starting at first position in the
source sequence would have ``start == 0`` and ``size == 1``.

As we can see on this example, ``start + size`` will give one more than the
zero-based end position. We can therefore manipulate ``start`` and
``start + size`` as python list slice boundaries.

For an inclusive end coordinate, we need to use ``end = start + size - 1``.
A 1-column wide alignment would have ``start == end``.
"""
MAFINDEX_VERSION = ...
class MafWriter(SequentialAlignmentWriter):
    """Accepts a MultipleSeqAlignment object, writes a MAF file."""
    def write_header(self): # -> None:
        """Write the MAF header."""
        ...
    
    def write_alignment(self, alignment):
        """Write a complete alignment to a MAF block.

        Writes every SeqRecord in a MultipleSeqAlignment object to its own
        MAF block (beginning with an 'a' line, containing 's' lines).
        """
        ...
    


def MafIterator(handle, seq_count=...):
    """Iterate over a MAF file handle as MultipleSeqAlignment objects.

    Iterates over lines in a MAF file-like object (handle), yielding
    MultipleSeqAlignment objects. SeqRecord IDs generally correspond to
    species names.
    """
    ...

class MafIndex:
    """Index for a MAF file.

    The index is a sqlite3 database that is built upon creation of the object
    if necessary, and queried when methods *search* or *get_spliced* are
    used.
    """
    def __init__(self, sqlite_file, maf_file, target_seqname) -> None:
        """Indexes or loads the index of a MAF file."""
        ...
    
    def search(self, starts, ends):
        """Search index database for MAF records overlapping ranges provided.

        Returns *MultipleSeqAlignment* results in order by start, then end, then
        internal offset field.

        *starts* should be a list of 0-based start coordinates of segments in the reference.
        *ends* should be the list of the corresponding segment ends
        (in the half-open UCSC convention:
        http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/).
        """
        ...
    
    def get_spliced(self, starts, ends, strand=...):
        """Return a multiple alignment of the exact sequence range provided.

        Accepts two lists of start and end positions on target_seqname, representing
        exons to be spliced in silico.  Returns a *MultipleSeqAlignment* of the
        desired sequences spliced together.

        *starts* should be a list of 0-based start coordinates of segments in the reference.
        *ends* should be the list of the corresponding segment ends
        (in the half-open UCSC convention:
        http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/).

        To ask for the alignment portion corresponding to the first 100
        nucleotides of the reference sequence, you would use
        ``search([0], [100])``
        """
        ...
    
    def __repr__(self): # -> str:
        """Return a string representation of the index."""
        ...
    
    def __len__(self): # -> int:
        """Return the number of records in the index."""
        ...
    


