"""
This type stub file was generated by pyright.
"""

from .Interfaces import SequenceIterator, SequenceWriter

"""Bio.SeqIO support for the "fasta" (aka FastA or Pearson) file format.

You are expected to use this module via the Bio.SeqIO functions.
"""
def SimpleFastaParser(handle): # -> Generator[tuple[Unknown, str], None, None]:
    """Iterate over Fasta records as string tuples.

    Arguments:
     - handle - input stream opened in text mode

    For each record a tuple of two strings is returned, the FASTA title
    line (without the leading '>' character), and the sequence (with any
    whitespace removed). The title line is not divided up into an
    identifier (the first word) and comment or description.

    >>> with open("Fasta/dups.fasta") as handle:
    ...     for values in SimpleFastaParser(handle):
    ...         print(values)
    ...
    ('alpha', 'ACGTA')
    ('beta', 'CGTC')
    ('gamma', 'CCGCC')
    ('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
    ('delta', 'CGCGC')

    """
    ...

def FastaTwoLineParser(handle): # -> Generator[tuple[Unbound | Unknown, Unknown], None, None]:
    """Iterate over no-wrapping Fasta records as string tuples.

    Arguments:
     - handle - input stream opened in text mode

    Functionally the same as SimpleFastaParser but with a strict
    interpretation of the FASTA format as exactly two lines per
    record, the greater-than-sign identifier with description,
    and the sequence with no line wrapping.

    Any line wrapping will raise an exception, as will excess blank
    lines (other than the special case of a zero-length sequence
    as the second line of a record).

    Examples
    --------
    This file uses two lines per FASTA record:

    >>> with open("Fasta/aster_no_wrap.pro") as handle:
    ...     for title, seq in FastaTwoLineParser(handle):
    ...         print("%s = %s..." % (title, seq[:3]))
    ...
    gi|3298468|dbj|BAA31520.1| SAMIPF = GGH...

    This equivalent file uses line wrapping:

    >>> with open("Fasta/aster.pro") as handle:
    ...     for title, seq in FastaTwoLineParser(handle):
    ...         print("%s = %s..." % (title, seq[:3]))
    ...
    Traceback (most recent call last):
       ...
    ValueError: Expected FASTA record starting with '>' character. Perhaps this file is using FASTA line wrapping? Got: 'MTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI'

    """
    ...

class FastaIterator(SequenceIterator):
    """Parser for Fasta files."""
    def __init__(self, source, alphabet=..., title2ids=...) -> None:
        """Iterate over Fasta records as SeqRecord objects.

        Arguments:
         - source - input stream opened in text mode, or a path to a file
         - alphabet - optional alphabet, not used. Leave as None.
         - title2ids - A function that, when given the title of the FASTA
           file (without the beginning >), will return the id, name and
           description (in that order) for the record as a tuple of strings.
           If this is not given, then the entire title line will be used
           as the description, and the first word as the id and name.

        By default this will act like calling Bio.SeqIO.parse(handle, "fasta")
        with no custom handling of the title lines:

        >>> with open("Fasta/dups.fasta") as handle:
        ...     for record in FastaIterator(handle):
        ...         print(record.id)
        ...
        alpha
        beta
        gamma
        alpha
        delta

        However, you can supply a title2ids function to alter this:

        >>> def take_upper(title):
        ...     return title.split(None, 1)[0].upper(), "", title
        >>> with open("Fasta/dups.fasta") as handle:
        ...     for record in FastaIterator(handle, title2ids=take_upper):
        ...         print(record.id)
        ...
        ALPHA
        BETA
        GAMMA
        ALPHA
        DELTA

        """
        ...
    
    def parse(self, handle): # -> Generator[SeqRecord, None, None]:
        """Start parsing the file, and return a SeqRecord generator."""
        ...
    
    def iterate(self, handle): # -> Generator[SeqRecord, None, None]:
        """Parse the file and generate SeqRecord objects."""
        ...
    


class FastaTwoLineIterator(SequenceIterator):
    """Parser for Fasta files with exactly two lines per record."""
    def __init__(self, source) -> None:
        """Iterate over two-line Fasta records (as SeqRecord objects).

        Arguments:
         - source - input stream opened in text mode, or a path to a file

        This uses a strict interpretation of the FASTA as requiring
        exactly two lines per record (no line wrapping).

        Only the default title to ID/name/description parsing offered
        by the relaxed FASTA parser is offered.
        """
        ...
    
    def parse(self, handle): # -> Generator[SeqRecord, None, None]:
        """Start parsing the file, and return a SeqRecord generator."""
        ...
    
    def iterate(self, handle): # -> Generator[SeqRecord, None, None]:
        """Parse the file and generate SeqRecord objects."""
        ...
    


class FastaWriter(SequenceWriter):
    """Class to write Fasta format files (OBSOLETE).

    Please use the ``as_fasta`` function instead, or the top level
    ``Bio.SeqIO.write()`` function instead using ``format="fasta"``.
    """
    def __init__(self, target, wrap=..., record2title=...) -> None:
        """Create a Fasta writer (OBSOLETE).

        Arguments:
         - target - Output stream opened in text mode, or a path to a file.
         - wrap -   Optional line length used to wrap sequence lines.
           Defaults to wrapping the sequence at 60 characters
           Use zero (or None) for no wrapping, giving a single
           long line for the sequence.
         - record2title - Optional function to return the text to be
           used for the title line of each record.  By default
           a combination of the record.id and record.description
           is used.  If the record.description starts with the
           record.id, then just the record.description is used.

        You can either use::

            handle = open(filename, "w")
            writer = FastaWriter(handle)
            writer.write_file(myRecords)
            handle.close()

        Or, follow the sequential file writer system, for example::

            handle = open(filename, "w")
            writer = FastaWriter(handle)
            writer.write_header() # does nothing for Fasta files
            ...
            Multiple writer.write_record() and/or writer.write_records() calls
            ...
            writer.write_footer() # does nothing for Fasta files
            handle.close()

        """
        ...
    
    def write_record(self, record):
        """Write a single Fasta record to the file."""
        ...
    


class FastaTwoLineWriter(FastaWriter):
    """Class to write 2-line per record Fasta format files (OBSOLETE).

    This means we write the sequence information  without line
    wrapping, and will always write a blank line for an empty
    sequence.

    Please use the ``as_fasta_2line`` function instead, or the top level
    ``Bio.SeqIO.write()`` function instead using ``format="fasta"``.
    """
    def __init__(self, handle, record2title=...) -> None:
        """Create a 2-line per record Fasta writer (OBSOLETE).

        Arguments:
         - handle - Handle to an output file, e.g. as returned
           by open(filename, "w")
         - record2title - Optional function to return the text to be
           used for the title line of each record.  By default
           a combination of the record.id and record.description
           is used.  If the record.description starts with the
           record.id, then just the record.description is used.

        You can either use::

            handle = open(filename, "w")
            writer = FastaWriter(handle)
            writer.write_file(myRecords)
            handle.close()

        Or, follow the sequential file writer system, for example::

            handle = open(filename, "w")
            writer = FastaWriter(handle)
            writer.write_header() # does nothing for Fasta files
            ...
            Multiple writer.write_record() and/or writer.write_records() calls
            ...
            writer.write_footer() # does nothing for Fasta files
            handle.close()

        """
        ...
    


def as_fasta(record): # -> str:
    """Turn a SeqRecord into a FASTA formatted string.

    This is used internally by the SeqRecord's .format("fasta")
    method and by the SeqIO.write(..., ..., "fasta") function.
    """
    ...

def as_fasta_2line(record): # -> str:
    """Turn a SeqRecord into a two-line FASTA formatted string.

    This is used internally by the SeqRecord's .format("fasta-2line")
    method and by the SeqIO.write(..., ..., "fasta-2line") function.
    """
    ...

if __name__ == "__main__":
    ...
