"""
This type stub file was generated by pyright.
"""

"""Parser for ACE files output by PHRAP.

Written by Frank Kauff (fkauff@duke.edu) and
Cymon J. Cox (cymon@duke.edu)

Usage:

There are two ways of reading an ace file:

1. The function 'read' reads the whole file at once;
2. The function 'parse' reads the file contig after contig.

First option, parse whole ace file at once::

        from Bio.Sequencing import Ace
        acefilerecord = Ace.read(open('my_ace_file.ace'))

This gives you:
 - acefilerecord.ncontigs (the number of contigs in the ace file)
 - acefilerecord.nreads (the number of reads in the ace file)
 - acefilerecord.contigs[] (one instance of the Contig class for each contig)

The Contig class holds the info of the CO tag, CT and WA tags, and all the reads used
for this contig in a list of instances of the Read class, e.g.::

        contig3 = acefilerecord.contigs[2]
        read4 = contig3.reads[3]
        RD_of_read4 = read4.rd
        DS_of_read4 = read4.ds

CT, WA, RT tags from the end of the file can appear anywhere are automatically
sorted into the right place.

see _RecordConsumer for details.

The second option is to  iterate over the contigs of an ace file one by one
in the ususal way::

    from Bio.Sequencing import Ace
    contigs = Ace.parse(open('my_ace_file.ace'))
    for contig in contigs:
        print(contig.name)
        ...

Please note that for memory efficiency, when using the iterator approach, only one
contig is kept in memory at once.  However, there can be a footer to the ACE file
containing WA, CT, RT or WR tags which contain additional meta-data on the contigs.
Because the parser doesn't see this data until the final record, it cannot be added to
the appropriate records.  Instead these tags will be returned with the last contig record.
Thus an ace file does not entirerly suit the concept of iterating. If WA, CT, RT, WR tags
are needed, the 'read' function rather than the 'parse' function might be more appropriate.
"""
class rd:
    """RD (reads), store a read with its name, sequence etc.

    The location and strand each read is mapped to is held in the AF lines.
    """
    def __init__(self) -> None:
        """Initialize the class."""
        ...
    


class qa:
    """QA (read quality), including which part if any was used as the consensus."""
    def __init__(self, line=...) -> None:
        """Initialize the class."""
        ...
    


class ds:
    """DS lines, include file name of a read's chromatogram file."""
    def __init__(self, line=...) -> None:
        """Initialize the class."""
        ...
    


class af:
    """AF lines, define the location of the read within the contig.

    Note attribute coru is short for complemented (C) or uncomplemented (U),
    since the strand information is stored in an ACE file using either the
    C or U character.
    """
    def __init__(self, line=...) -> None:
        """Initialize the class."""
        ...
    


class bs:
    """BS (base segment), which read was chosen as the consensus at each position."""
    def __init__(self, line=...) -> None:
        """Initialize the class."""
        ...
    


class rt:
    """RT (transient read tags), generated by crossmatch and phrap."""
    def __init__(self, line=...) -> None:
        """Initialize the class."""
        ...
    


class ct:
    """CT (consensus tags)."""
    def __init__(self, line=...) -> None:
        """Initialize the class."""
        ...
    


class wa:
    """WA (whole assembly tag), holds the assembly program name, version, etc."""
    def __init__(self, line=...) -> None:
        """Initialize the class."""
        ...
    


class wr:
    """WR lines."""
    def __init__(self, line=...) -> None:
        """Initialize the class."""
        ...
    


class Reads:
    """Holds information about a read supporting an ACE contig."""
    def __init__(self, line=...) -> None:
        """Initialize the class."""
        ...
    


class Contig:
    """Holds information about a contig from an ACE record."""
    def __init__(self, line=...) -> None:
        """Initialize the class."""
        ...
    


def parse(source): # -> Generator[Contig, None, None]:
    """Iterate of ACE file contig by contig.

    Argument source is a file-like object or a path to a file.

    This function returns an iterator that allows you to iterate
    over the ACE file record by record::

        records = parse(source)
        for record in records:
            # do something with the record

    where each record is a Contig object.
    """
    ...

class ACEFileRecord:
    """Holds data of an ACE file."""
    def __init__(self) -> None:
        """Initialize the class."""
        ...
    
    def sort(self): # -> None:
        """Sorts wr, rt and ct tags into the appropriate contig / read instance, if possible."""
        ...
    


def read(handle): # -> ACEFileRecord:
    """Parse a full ACE file into a list of contigs."""
    ...

