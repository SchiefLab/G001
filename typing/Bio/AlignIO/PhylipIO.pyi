"""
This type stub file was generated by pyright.
"""

from .Interfaces import AlignmentIterator, SequentialAlignmentWriter

"""AlignIO support for "phylip" format from Joe Felsenstein's PHYLIP tools.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).

Support for "relaxed phylip" format is also provided. Relaxed phylip differs
from standard phylip format in the following ways:

 - No whitespace is allowed in the sequence ID.
 - No truncation is performed. Instead, sequence IDs are padded to the longest
   ID length, rather than 10 characters. A space separates the sequence
   identifier from the sequence.

Relaxed phylip is supported by RAxML and PHYML.

Note
====

In TREE_PUZZLE (Schmidt et al. 2003) and PHYML (Guindon and Gascuel 2003)
a dot/period (".") in a sequence is interpreted as meaning the same
character as in the first sequence.  The PHYLIP documentation from 3.3 to 3.69
http://evolution.genetics.washington.edu/phylip/doc/sequence.html says:

"a period was also previously allowed but it is no longer allowed,
because it sometimes is used in different senses in other programs"

Biopython 1.58 or later treats dots/periods in the sequence as invalid, both
for reading and writing. Older versions did nothing special with a dot/period.
"""
_PHYLIP_ID_WIDTH = ...
_NO_DOTS = ...
class PhylipWriter(SequentialAlignmentWriter):
    """Phylip alignment writer."""
    def write_alignment(self, alignment, id_width=...): # -> None:
        """Use this to write (another) single alignment to an open file.

        This code will write interlaced alignments (when the sequences are
        longer than 50 characters).

        Note that record identifiers are strictly truncated to id_width,
        defaulting to the value required to comply with the PHYLIP standard.

        For more information on the file format, please see:
        http://evolution.genetics.washington.edu/phylip/doc/sequence.html
        http://evolution.genetics.washington.edu/phylip/doc/main.html#inputfiles
        """
        ...
    


class PhylipIterator(AlignmentIterator):
    """Reads a Phylip alignment file returning a MultipleSeqAlignment iterator.

    Record identifiers are limited to at most 10 characters.

    It only copes with interlaced phylip files!  Sequential files won't work
    where the sequences are split over multiple lines.

    For more information on the file format, please see:
    http://evolution.genetics.washington.edu/phylip/doc/sequence.html
    http://evolution.genetics.washington.edu/phylip/doc/main.html#inputfiles
    """
    id_width = ...
    _header = ...
    def __next__(self): # -> MultipleSeqAlignment:
        """Parse the next alignment from the handle."""
        ...
    


class RelaxedPhylipWriter(PhylipWriter):
    """Relaxed Phylip format writer."""
    def write_alignment(self, alignment): # -> None:
        """Write a relaxed phylip alignment."""
        ...
    


class RelaxedPhylipIterator(PhylipIterator):
    """Relaxed Phylip format Iterator."""
    ...


class SequentialPhylipWriter(SequentialAlignmentWriter):
    """Sequential Phylip format Writer."""
    def write_alignment(self, alignment, id_width=...): # -> None:
        """Write a Phylip alignment to the handle."""
        ...
    


class SequentialPhylipIterator(PhylipIterator):
    """Sequential Phylip format Iterator.

    The sequential format carries the same restrictions as the normal
    interleaved one, with the difference being that the sequences are listed
    sequentially, each sequence written in its entirety before the start of
    the next. According to the PHYLIP documentation for input file
    formatting, newlines and spaces may optionally be entered at any point
    in the sequences.
    """
    _header = ...
    def __next__(self): # -> MultipleSeqAlignment:
        """Parse the next alignment from the handle."""
        ...
    


def sanitize_name(name, width=...):
    """Sanitise sequence identifier for output.

    Removes the banned characters "[]()" and replaces the characters ":;"
    with "|". The name is truncated to "width" characters if specified.
    """
    ...

