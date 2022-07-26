"""
This type stub file was generated by pyright.
"""

from .Interfaces import AlignmentIterator, SequentialAlignmentWriter

"""Bio.AlignIO support for "xmfa" output from Mauve/ProgressiveMauve.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).

For example, consider a progressiveMauve alignment file containing the following::

    #FormatVersion Mauve1
    #Sequence1File	a.fa
    #Sequence1Entry	1
    #Sequence1Format	FastA
    #Sequence2File	b.fa
    #Sequence2Entry	2
    #Sequence2Format	FastA
    #Sequence3File	c.fa
    #Sequence3Entry	3
    #Sequence3Format	FastA
    #BackboneFile	three.xmfa.bbcols
    > 1:0-0 + a.fa
    --------------------------------------------------------------------------------
    --------------------------------------------------------------------------------
    --------------------------------------------------------------------------------
    > 2:5417-5968 + b.fa
    TTTAAACATCCCTCGGCCCGTCGCCCTTTTATAATAGCAGTACGTGAGAGGAGCGCCCTAAGCTTTGGGAAATTCAAGC-
    --------------------------------------------------------------------------------
    CTGGAACGTACTTGCTGGTTTCGCTACTATTTCAAACAAGTTAGAGGCCGTTACCTCGGGCGAACGTATAAACCATTCTG
    > 3:9476-10076 - c.fa
    TTTAAACACCTTTTTGGATG--GCCCAGTTCGTTCAGTTGTG-GGGAGGAGATCGCCCCAAACGTATGGTGAGTCGGGCG
    TTTCCTATAGCTATAGGACCAATCCACTTACCATACGCCCGGCGTCGCCCAGTCCGGTTCGGTACCCTCCATGACCCACG
    ---------------------------------------------------------AAATGAGGGCCCAGGGTATGCTT
    =
    > 2:5969-6015 + b.fa
    -----------------------
    GGGCGAACGTATAAACCATTCTG
    > 3:9429-9476 - c.fa
    TTCGGTACCCTCCATGACCCACG
    AAATGAGGGCCCAGGGTATGCTT

This is a multiple sequence alignment with multiple aligned sections, so you
would probably load this using the Bio.AlignIO.parse() function:

    >>> from Bio import AlignIO
    >>> align = AlignIO.parse("Mauve/simple_short.xmfa", "mauve")
    >>> alignments = list(align)
    >>> for aln in alignments:
    ...     print(aln)
    ...
    Alignment with 3 rows and 240 columns
    --------------------------------------------...--- a.fa
    TTTAAACATCCCTCGGCCCGTCGCCCTTTTATAATAGCAGTACG...CTG b.fa/5416-5968
    TTTAAACACCTTTTTGGATG--GCCCAGTTCGTTCAGTTGTG-G...CTT c.fa/9475-10076
    Alignment with 2 rows and 46 columns
    -----------------------GGGCGAACGTATAAACCATTCTG b.fa/5968-6015
    TTCGGTACCCTCCATGACCCACGAAATGAGGGCCCAGGGTATGCTT c.fa/9428-9476

Additional information is extracted from the XMFA file and available through
the annotation attribute of each record::

    >>> for record in alignments[0]:
    ...     print(record.id, len(record))
    ...     print("  start: %d, end: %d, strand: %d" %(
    ...         record.annotations['start'], record.annotations['end'],
    ...         record.annotations['strand']))
    ...
    a.fa 240
      start: 0, end: 0, strand: 1
    b.fa/5416-5968 240
      start: 5416, end: 5968, strand: 1
    c.fa/9475-10076 240
      start: 9475, end: 10076, strand: -1

"""
XMFA_HEADER_REGEX = ...
XMFA_HEADER_REGEX_BIOPYTHON = ...
ID_LINE_FMT = ...
class MauveWriter(SequentialAlignmentWriter):
    """Mauve/XMFA alignment writer."""
    def __init__(self, *args, **kwargs) -> None:
        """Initialize the class."""
        ...
    
    def write_alignment(self, alignment): # -> None:
        """Use this to write (another) single alignment to an open file.

        Note that sequences and their annotation are recorded
        together (rather than having a block of annotation followed
        by a block of aligned sequences).
        """
        ...
    


class MauveIterator(AlignmentIterator):
    """Mauve xmfa alignment iterator."""
    _ids = ...
    def __next__(self): # -> MultipleSeqAlignment:
        """Parse the next alignment from the handle."""
        ...
    


