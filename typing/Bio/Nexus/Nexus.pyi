"""
This type stub file was generated by pyright.
"""

"""Nexus class. Parse the contents of a NEXUS file.

Based upon 'NEXUS: An extensible file format for systematic information'
Maddison, Swofford, Maddison. 1997. Syst. Biol. 46(4):590-621
"""
INTERLEAVE = ...
SPECIAL_COMMANDS = ...
KNOWN_NEXUS_BLOCKS = ...
PUNCTUATION = ...
MRBAYESSAFE = ...
WHITESPACE = ...
SPECIALCOMMENTS = ...
CHARSET = ...
TAXSET = ...
CODONPOSITIONS = ...
DEFAULTNEXUS = ...
class NexusError(Exception):
    """Provision for the management of Nexus exceptions."""
    ...


class CharBuffer:
    """Helps reading NEXUS-words and characters from a buffer (semi-PRIVATE).

    This class is not intended for public use (any more).
    """
    def __init__(self, string) -> None:
        """Initialize the class."""
        ...
    
    def peek(self): # -> None:
        """Return the first character from the buffer."""
        ...
    
    def peek_nonwhitespace(self): # -> str | None:
        """Return the first character from the buffer, do not include spaces."""
        ...
    
    def __next__(self): # -> None:
        """Iterate over NEXUS characters in the file."""
        ...
    
    def next_nonwhitespace(self): # -> None:
        """Check for next non whitespace character in NEXUS file."""
        ...
    
    def skip_whitespace(self): # -> None:
        """Skip whitespace characters in NEXUS file."""
        ...
    
    def next_until(self, target): # -> str | None:
        """Iterate over the NEXUS file until a target character is reached."""
        ...
    
    def peek_word(self, word):
        """Return a word stored in the buffer."""
        ...
    
    def next_word(self): # -> str | None:
        """Return the next NEXUS word from a string.

        This deals with single and double quotes, whitespace and punctuation.
        """
        ...
    
    def rest(self): # -> str:
        """Return the rest of the string without parsing."""
        ...
    


class StepMatrix:
    """Calculate a stepmatrix for weighted parsimony.

    See :
    COMBINATORIAL WEIGHTS IN PHYLOGENETIC ANALYSIS - A STATISTICAL PARSIMONY PROCEDURE
    Wheeler (1990), Cladistics 6:269-275.
    """
    def __init__(self, symbols, gap) -> None:
        """Initialize the class."""
        ...
    
    def set(self, x, y, value): # -> None:
        """Set a given value in the matrix's position."""
        ...
    
    def add(self, x, y, value): # -> None:
        """Add the given value to existing, in matrix's position."""
        ...
    
    def sum(self):
        """Calculate the associations, makes matrix of associations."""
        ...
    
    def transformation(self): # -> Self@StepMatrix:
        """Calculate the transformation matrix.

        Normalizes the columns of the matrix of associations.
        """
        ...
    
    def weighting(self): # -> Self@StepMatrix:
        """Calculate the Phylogenetic weight matrix.

        Constructed from the logarithmic transformation of the
        transformation matrix.
        """
        ...
    
    def smprint(self, name=...): # -> str:
        """Print a stepmatrix."""
        ...
    


def safename(name, mrbayes=...): # -> str:
    """Return a taxon identifier according to NEXUS standard.

    Wrap quotes around names with punctuation or whitespace, and double
    single quotes.

    mrbayes=True: write names without quotes, whitespace or punctuation
    for the mrbayes software package.
    """
    ...

def quotestrip(word): # -> None:
    """Remove quotes and/or double quotes around identifiers."""
    ...

def get_start_end(sequence, skiplist=...): # -> tuple[None, None] | tuple[Literal[-1], Literal[-1]] | tuple[int, int]:
    """Return position of first and last character which is not in skiplist.

    Skiplist defaults to ['-','?'].
    """
    ...

def combine(matrices): # -> None:
    """Combine matrices in [(name,nexus-instance),...] and return new nexus instance.

    combined_matrix=combine([(name1,nexus_instance1),(name2,nexus_instance2),...]
    Character sets, character partitions and taxon sets are prefixed, readjusted
    and present in the combined matrix.
    """
    ...

class Commandline:
    """Represent a commandline as command and options."""
    def __init__(self, line, title) -> None:
        """Initialize the class."""
        ...
    


class Block:
    """Represent a NEXUS block with block name and list of commandlines."""
    def __init__(self, title=...) -> None:
        """Initialize the class."""
        ...
    


class Nexus:
    """Create the Nexus class, main class for the management of Nexus files."""
    def __init__(self, input=...) -> None:
        """Initialize the class."""
        ...
    
    def get_original_taxon_order(self): # -> list[Unknown]:
        """Included for backwards compatibility (DEPRECATED)."""
        ...
    
    def set_original_taxon_order(self, value): # -> None:
        """Included for backwards compatibility (DEPRECATED)."""
        ...
    
    original_taxon_order = ...
    def read(self, input):
        """Read and parse NEXUS input (a filename, file-handle, or string)."""
        ...
    
    def write_nexus_data_partitions(self, matrix=..., filename=..., blocksize=..., interleave=..., exclude=..., delete=..., charpartition=..., comment=..., mrbayes=...):
        """Write a nexus file for each partition in charpartition.

        Only non-excluded characters and non-deleted taxa are included,
        just the data block is written.
        """
        ...
    
    def write_nexus_data(self, filename=..., matrix=..., exclude=..., delete=..., blocksize=..., interleave=..., interleave_by_partition=..., comment=..., omit_NEXUS=..., append_sets=..., mrbayes=..., codons_block=...):
        """Write a nexus file with data and sets block to a file or handle.

        Character sets and partitions are appended by default, and are
        adjusted according to excluded characters (i.e. character sets
        still point to the same sites (not necessarily same positions),
        without including the deleted characters.

        - filename - Either a filename as a string (which will be opened,
          written to and closed), or a handle object (which will
          be written to but NOT closed).
        - interleave_by_partition - Optional name of partition (string)
        - omit_NEXUS - Boolean.  If true, the '#NEXUS' line normally at the
          start of the file is omitted.

        Returns the filename/handle used to write the data.
        """
        ...
    
    def append_sets(self, exclude=..., delete=..., mrbayes=..., include_codons=..., codons_only=...):
        """Return a sets block."""
        ...
    
    def export_fasta(self, filename=..., width=...): # -> str | Any:
        """Write matrix into a fasta file."""
        ...
    
    def export_phylip(self, filename=...): # -> str | Any:
        """Write matrix into a PHYLIP file.

        Note that this writes a relaxed PHYLIP format file, where the names
        are not truncated, nor checked for invalid characters.
        """
        ...
    
    def constant(self, matrix=..., delete=..., exclude=...):
        """Return a list with all constant characters."""
        ...
    
    def cstatus(self, site, delete=..., narrow=...):
        """Summarize character.

        narrow=True:  paup-mode (a c ? --> ac; ? ? ? --> ?)
        narrow=false:           (a c ? --> a c g t -; ? ? ? --> a c g t -)
        """
        ...
    
    def weighted_stepmatrix(self, name=..., exclude=..., delete=...): # -> str:
        """Calculate a stepmatrix for weighted parsimony.

        See Wheeler (1990), Cladistics 6:269-275 and
        Felsenstein (1981), Biol. J. Linn. Soc. 16:183-196
        """
        ...
    
    def crop_matrix(self, matrix=..., delete=..., exclude=...):
        """Return a matrix without deleted taxa and excluded characters."""
        ...
    
    def bootstrap(self, matrix=..., delete=..., exclude=...):
        """Return a bootstrapped matrix."""
        ...
    
    def add_sequence(self, name, sequence): # -> None:
        """Add a sequence (string) to the matrix."""
        ...
    
    def insert_gap(self, pos, n=..., leftgreedy=...): # -> None:
        """Add a gap into the matrix and adjust charsets and partitions.

        pos=0: first position
        pos=nchar: last position
        """
        ...
    
    def invert(self, charlist): # -> list[int]:
        """Return all character indices that are not in charlist."""
        ...
    
    def gaponly(self, include_missing=...): # -> list[int]:
        """Return gap-only sites."""
        ...
    
    def terminal_gap_to_missing(self, missing=..., skip_n=...): # -> None:
        """Replace all terminal gaps with missing character.

        Mixtures like ???------??------- are properly resolved.
        """
        ...
    


if __name__ == "__main__":
    ...
