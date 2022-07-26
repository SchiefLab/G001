"""
This type stub file was generated by pyright.
"""

import contextlib
import collections.abc
from abc import ABC, abstractmethod

"""Code for more fancy file handles.

Bio.File defines private classes used in Bio.SeqIO and Bio.SearchIO for
indexing files. These are not intended for direct use.
"""
@contextlib.contextmanager
def as_handle(handleish, mode=..., **kwargs): # -> Generator[IO[Any] | Unknown, None, None]:
    r"""Context manager to ensure we are using a handle.

    Context manager for arguments that can be passed to SeqIO and AlignIO read, write,
    and parse methods: either file objects or path-like objects (strings, pathlib.Path
    instances, or more generally, anything that can be handled by the builtin 'open'
    function).

    When given a path-like object, returns an open file handle to that path, with provided
    mode, which will be closed when the manager exits.

    All other inputs are returned, and are *not* closed.

    Arguments:
     - handleish  - Either a file handle or path-like object (anything which can be
                    passed to the builtin 'open' function, such as str, bytes,
                    pathlib.Path, and os.DirEntry objects)
     - mode       - Mode to open handleish (used only if handleish is a string)
     - kwargs     - Further arguments to pass to open(...)

    Examples
    --------
    >>> from Bio import File
    >>> import os
    >>> with File.as_handle('seqs.fasta', 'w') as fp:
    ...     fp.write('>test\nACGT')
    ...
    10
    >>> fp.closed
    True

    >>> handle = open('seqs.fasta', 'w')
    >>> with File.as_handle(handle) as fp:
    ...     fp.write('>test\nACGT')
    ...
    10
    >>> fp.closed
    False
    >>> fp.close()
    >>> os.remove("seqs.fasta")  # tidy up

    """
    ...

class _IndexedSeqFileProxy(ABC):
    """Abstract base class for file format specific random access (PRIVATE).

    This is subclasses in both Bio.SeqIO for indexing as SeqRecord
    objects, and in Bio.SearchIO for indexing QueryResult objects.

    Subclasses for each file format should define '__iter__', 'get'
    and optionally 'get_raw' methods.
    """
    @abstractmethod
    def __iter__(self):
        """Return (identifier, offset, length in bytes) tuples.

        The length can be zero where it is not implemented or not
        possible for a particular file format.
        """
        ...
    
    @abstractmethod
    def get(self, offset):
        """Return parsed object for this entry."""
        ...
    
    def get_raw(self, offset):
        """Return the raw record from the file as a bytes string (if implemented).

        If the key is not found, a KeyError exception is raised.

        This may not have been implemented for all file formats.
        """
        ...
    


class _IndexedSeqFileDict(collections.abc.Mapping):
    """Read only dictionary interface to a sequential record file.

    This code is used in both Bio.SeqIO for indexing as SeqRecord
    objects, and in Bio.SearchIO for indexing QueryResult objects.

    Keeps the keys and associated file offsets in memory, reads the file
    to access entries as objects parsing them on demand. This approach
    is memory limited, but will work even with millions of records.

    Note duplicate keys are not allowed. If this happens, a ValueError
    exception is raised.

    As used in Bio.SeqIO, by default the SeqRecord's id string is used
    as the dictionary key. In Bio.SearchIO, the query's id string is
    used. This can be changed by supplying an optional key_function,
    a callback function which will be given the record id and must
    return the desired key. For example, this allows you to parse
    NCBI style FASTA identifiers, and extract the GI number to use
    as the dictionary key.

    Note that this dictionary is essentially read only. You cannot
    add or change values, pop values, nor clear the dictionary.
    """
    def __init__(self, random_access_proxy, key_function, repr, obj_repr) -> None:
        """Initialize the class."""
        ...
    
    def __repr__(self): # -> Unknown:
        """Return a string representation of the File object."""
        ...
    
    def __str__(self) -> str:
        """Create a string representation of the File object."""
        ...
    
    def __len__(self): # -> int:
        """Return the number of records."""
        ...
    
    def __iter__(self): # -> Iterator[Unknown]:
        """Iterate over the keys."""
        ...
    
    def __getitem__(self, key):
        """Return record for the specified key."""
        ...
    
    def get_raw(self, key):
        """Return the raw record from the file as a bytes string.

        If the key is not found, a KeyError exception is raised.
        """
        ...
    
    def close(self): # -> None:
        """Close the file handle being used to read the data.

        Once called, further use of the index won't work. The sole purpose
        of this method is to allow explicit handle closure - for example
        if you wish to delete the file, on Windows you must first close
        all open handles to that file.
        """
        ...
    


class _SQLiteManySeqFilesDict(_IndexedSeqFileDict):
    """Read only dictionary interface to many sequential record files.

    This code is used in both Bio.SeqIO for indexing as SeqRecord
    objects, and in Bio.SearchIO for indexing QueryResult objects.

    Keeps the keys, file-numbers and offsets in an SQLite database. To access
    a record by key, reads from the offset in the appropriate file and then
    parses the record into an object.

    There are OS limits on the number of files that can be open at once,
    so a pool are kept. If a record is required from a closed file, then
    one of the open handles is closed first.
    """
    def __init__(self, index_filename, filenames, proxy_factory, fmt, key_function, repr, max_open=...) -> None:
        """Initialize the class."""
        ...
    
    def __repr__(self):
        ...
    
    def __contains__(self, key): # -> bool:
        ...
    
    def __len__(self): # -> int:
        """Return the number of records indexed."""
        ...
    
    def __iter__(self): # -> Generator[str, None, None]:
        """Iterate over the keys."""
        ...
    
    def __getitem__(self, key):
        """Return record for the specified key."""
        ...
    
    def get_raw(self, key):
        """Return the raw record from the file as a bytes string.

        If the key is not found, a KeyError exception is raised.
        """
        ...
    
    def close(self): # -> None:
        """Close any open file handles."""
        ...
    


