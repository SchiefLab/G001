"""
This type stub file was generated by pyright.
"""

"""Hold GenBank data in a straightforward format.

Classes:
 - Record - All of the information in a GenBank record.
 - Reference - hold reference data for a record.
 - Feature - Hold the information in a Feature Table.
 - Qualifier - Qualifiers on a Feature.

"""
class Record:
    """Hold GenBank information in a format similar to the original record.

    The Record class is meant to make data easy to get to when you are
    just interested in looking at GenBank data.

    Attributes:
     - locus - The name specified after the LOCUS keyword in the GenBank
       record. This may be the accession number, or a clone id or something else.
     - size - The size of the record.
     - residue_type - The type of residues making up the sequence in this
       record. Normally something like RNA, DNA or PROTEIN, but may be as
       esoteric as 'ss-RNA circular'.
     - data_file_division - The division this record is stored under in
       GenBank (ie. PLN -> plants; PRI -> humans, primates; BCT -> bacteria...)
     - date - The date of submission of the record, in a form like '28-JUL-1998'
     - accession - list of all accession numbers for the sequence.
     - nid - Nucleotide identifier number.
     - pid - Proteint identifier number
     - version - The accession number + version (ie. AB01234.2)
     - db_source - Information about the database the record came from
     - gi - The NCBI gi identifier for the record.
     - keywords - A list of keywords related to the record.
     - segment - If the record is one of a series, this is info about which
       segment this record is (something like '1 of 6').
     - source - The source of material where the sequence came from.
     - organism - The genus and species of the organism (ie. 'Homo sapiens')
     - taxonomy - A listing of the taxonomic classification of the organism,
       starting general and getting more specific.
     - references - A list of Reference objects.
     - comment - Text with any kind of comment about the record.
     - features - A listing of Features making up the feature table.
     - base_counts - A string with the counts of bases for the sequence.
     - origin - A string specifying info about the origin of the sequence.
     - sequence - A string with the sequence itself.
     - contig - A string of location information for a CONTIG in a RefSeq file
     - project - The genome sequencing project numbers
       (will be replaced by the dblink cross-references in 2009).
     - dblinks - The genome sequencing project number(s) and other links.
       (will replace the project information in 2009).

    """
    GB_LINE_LENGTH = ...
    GB_BASE_INDENT = ...
    GB_FEATURE_INDENT = ...
    GB_INTERNAL_INDENT = ...
    GB_OTHER_INTERNAL_INDENT = ...
    GB_FEATURE_INTERNAL_INDENT = ...
    GB_SEQUENCE_INDENT = ...
    BASE_FORMAT = ...
    INTERNAL_FORMAT = ...
    OTHER_INTERNAL_FORMAT = ...
    BASE_FEATURE_FORMAT = ...
    INTERNAL_FEATURE_FORMAT = ...
    SEQUENCE_FORMAT = ...
    def __init__(self) -> None:
        """Initialize the class."""
        ...
    
    def __str__(self) -> str:
        """Provide a GenBank formatted output option for a Record.

        The objective of this is to provide an easy way to read in a GenBank
        record, modify it somehow, and then output it in 'GenBank format.'
        We are striving to make this work so that a parsed Record that is
        output using this function will look exactly like the original
        record.

        Much of the output is based on format description info at:

        ftp://ncbi.nlm.nih.gov/genbank/gbrel.txt
        """
        ...
    


class Reference:
    """Hold information from a GenBank reference.

    Attributes:
     - number - The number of the reference in the listing of references.
     - bases - The bases in the sequence the reference refers to.
     - authors - String with all of the authors.
     - consrtm - Consortium the authors belong to.
     - title - The title of the reference.
     - journal - Information about the journal where the reference appeared.
     - medline_id - The medline id for the reference.
     - pubmed_id - The pubmed_id for the reference.
     - remark - Free-form remarks about the reference.

    """
    def __init__(self) -> None:
        """Initialize the class."""
        ...
    
    def __str__(self) -> str:
        """Convert the reference to a GenBank format string."""
        ...
    


class Feature:
    """Hold information about a Feature in the Feature Table of GenBank record.

    Attributes:
     - key - The key name of the featue (ie. source)
     - location - The string specifying the location of the feature.
     - qualfiers - A list of Qualifier objects in the feature.

    """
    def __init__(self, key=..., location=...) -> None:
        """Initialize the class."""
        ...
    
    def __repr__(self): # -> str:
        """Representation of the object for debugging or logging."""
        ...
    
    def __str__(self) -> str:
        """Return feature as a GenBank format string."""
        ...
    


class Qualifier:
    """Hold information about a qualifier in a GenBank feature.

    Attributes:
     - key - The key name of the qualifier (ie. /organism=)
     - value - The value of the qualifier ("Dictyostelium discoideum").

    """
    def __init__(self, key=..., value=...) -> None:
        """Initialize the class."""
        ...
    
    def __repr__(self): # -> str:
        """Representation of the object for debugging or logging."""
        ...
    
    def __str__(self) -> str:
        """Return feature qualifier as a GenBank format string."""
        ...
    


