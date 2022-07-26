"""
This type stub file was generated by pyright.
"""

from .Interfaces import SequenceIterator, SequenceWriter

"""Bio.SeqIO support for the "genbank" and "embl" file formats.

You are expected to use this module via the Bio.SeqIO functions.
Note that internally this module calls Bio.GenBank to do the actual
parsing of GenBank, EMBL and IMGT files.

See Also:
International Nucleotide Sequence Database Collaboration
http://www.insdc.org/

GenBank
http://www.ncbi.nlm.nih.gov/Genbank/

EMBL Nucleotide Sequence Database
http://www.ebi.ac.uk/embl/

DDBJ (DNA Data Bank of Japan)
http://www.ddbj.nig.ac.jp/

IMGT (use a variant of EMBL format with longer feature indents)
http://imgt.cines.fr/download/LIGM-DB/userman_doc.html
http://imgt.cines.fr/download/LIGM-DB/ftable_doc.html
http://www.ebi.ac.uk/imgt/hla/docs/manual.html

"""
class GenBankIterator(SequenceIterator):
    """Parser for GenBank files."""
    def __init__(self, source) -> None:
        """Break up a Genbank file into SeqRecord objects.

        Argument source is a file-like object opened in text mode or a path to a file.
        Every section from the LOCUS line to the terminating // becomes
        a single SeqRecord with associated annotation and features.

        Note that for genomes or chromosomes, there is typically only
        one record.

        This gets called internally by Bio.SeqIO for the GenBank file format:

        >>> from Bio import SeqIO
        >>> for record in SeqIO.parse("GenBank/cor6_6.gb", "gb"):
        ...     print(record.id)
        ...
        X55053.1
        X62281.1
        M81224.1
        AJ237582.1
        L31939.1
        AF297471.1

        Equivalently,

        >>> with open("GenBank/cor6_6.gb") as handle:
        ...     for record in GenBankIterator(handle):
        ...         print(record.id)
        ...
        X55053.1
        X62281.1
        M81224.1
        AJ237582.1
        L31939.1
        AF297471.1

        """
        ...
    
    def parse(self, handle): # -> Generator[SeqRecord, None, None]:
        """Start parsing the file, and return a SeqRecord generator."""
        ...
    


class EmblIterator(SequenceIterator):
    """Parser for EMBL files."""
    def __init__(self, source) -> None:
        """Break up an EMBL file into SeqRecord objects.

        Argument source is a file-like object opened in text mode or a path to a file.
        Every section from the LOCUS line to the terminating // becomes
        a single SeqRecord with associated annotation and features.

        Note that for genomes or chromosomes, there is typically only
        one record.

        This gets called internally by Bio.SeqIO for the EMBL file format:

        >>> from Bio import SeqIO
        >>> for record in SeqIO.parse("EMBL/epo_prt_selection.embl", "embl"):
        ...     print(record.id)
        ...
        A00022.1
        A00028.1
        A00031.1
        A00034.1
        A00060.1
        A00071.1
        A00072.1
        A00078.1
        CQ797900.1

        Equivalently,

        >>> with open("EMBL/epo_prt_selection.embl") as handle:
        ...     for record in EmblIterator(handle):
        ...         print(record.id)
        ...
        A00022.1
        A00028.1
        A00031.1
        A00034.1
        A00060.1
        A00071.1
        A00072.1
        A00078.1
        CQ797900.1

        """
        ...
    
    def parse(self, handle): # -> Generator[SeqRecord, None, None]:
        """Start parsing the file, and return a SeqRecord generator."""
        ...
    


class ImgtIterator(SequenceIterator):
    """Parser for IMGT files."""
    def __init__(self, source) -> None:
        """Break up an IMGT file into SeqRecord objects.

        Argument source is a file-like object opened in text mode or a path to a file.
        Every section from the LOCUS line to the terminating // becomes
        a single SeqRecord with associated annotation and features.

        Note that for genomes or chromosomes, there is typically only
        one record.
        """
        ...
    
    def parse(self, handle): # -> Generator[SeqRecord, None, None]:
        """Start parsing the file, and return a SeqRecord generator."""
        ...
    


class GenBankCdsFeatureIterator(SequenceIterator):
    """Parser for GenBank files, creating a SeqRecord for each CDS feature."""
    def __init__(self, source) -> None:
        """Break up a Genbank file into SeqRecord objects for each CDS feature.

        Argument source is a file-like object opened in text mode or a path to a file.

        Every section from the LOCUS line to the terminating // can contain
        many CDS features.  These are returned as with the stated amino acid
        translation sequence (if given).
        """
        ...
    
    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        ...
    


class EmblCdsFeatureIterator(SequenceIterator):
    """Parser for EMBL files, creating a SeqRecord for each CDS feature."""
    def __init__(self, source) -> None:
        """Break up a EMBL file into SeqRecord objects for each CDS feature.

        Argument source is a file-like object opened in text mode or a path to a file.

        Every section from the LOCUS line to the terminating // can contain
        many CDS features.  These are returned as with the stated amino acid
        translation sequence (if given).
        """
        ...
    
    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        ...
    


class _InsdcWriter(SequenceWriter):
    """Base class for GenBank and EMBL writers (PRIVATE)."""
    MAX_WIDTH = ...
    QUALIFIER_INDENT = ...
    QUALIFIER_INDENT_STR = ...
    QUALIFIER_INDENT_TMP = ...
    FTQUAL_NO_QUOTE = ...


class GenBankWriter(_InsdcWriter):
    """GenBank writer."""
    HEADER_WIDTH = ...
    QUALIFIER_INDENT = ...
    STRUCTURED_COMMENT_START = ...
    STRUCTURED_COMMENT_END = ...
    STRUCTURED_COMMENT_DELIM = ...
    LETTERS_PER_LINE = ...
    SEQUENCE_INDENT = ...
    def write_record(self, record):
        """Write a single record to the output file."""
        ...
    


class EmblWriter(_InsdcWriter):
    """EMBL writer."""
    HEADER_WIDTH = ...
    QUALIFIER_INDENT = ...
    QUALIFIER_INDENT_STR = ...
    QUALIFIER_INDENT_TMP = ...
    FEATURE_HEADER = ...
    LETTERS_PER_BLOCK = ...
    BLOCKS_PER_LINE = ...
    LETTERS_PER_LINE = ...
    POSITION_PADDING = ...
    def write_record(self, record):
        """Write a single record to the output file."""
        ...
    


class ImgtWriter(EmblWriter):
    """IMGT writer (EMBL format variant)."""
    HEADER_WIDTH = ...
    QUALIFIER_INDENT = ...
    QUALIFIER_INDENT_STR = ...
    QUALIFIER_INDENT_TMP = ...
    FEATURE_HEADER = ...


if __name__ == "__main__":
    ...
