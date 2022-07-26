"""
This type stub file was generated by pyright.
"""

from Bio.PDB.Atom import Atom
from typing import Dict, List, Optional, TYPE_CHECKING, TextIO, Tuple, Union
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain

"""Classes to support internal coordinates for protein structures.

Internal coordinates comprise Psi, Phi and Omega dihedral angles along the
protein backbone, Chi angles along the sidechains, and all 3-atom angles and
bond lengths comprising a protein chain.  These routines can compute internal
coordinates from atom XYZ coordinates, and compute atom XYZ coordinates from
internal coordinates.

Internal coordinates are defined on sequences of atoms which span
residues or follow accepted nomenclature along sidechains.  To manage these
sequences and support Biopython's disorder mechanisms, AtomKey specifiers are
implemented to capture residue, atom and variant identification in a single
object.  A Hedron object is specified as three sequential AtomKeys, comprising
two bond lengths and the bond angle between them.  A Dihedron consists of four
sequential AtomKeys, linking two Hedra with a dihedral angle between them.

A Protein Internal Coordinate (.pic) file format is defined to capture
sufficient detail to reproduce a PDB file from chain starting coordinates
(first residue N, Ca, C XYZ coordinates) and remaining internal coordinates.
These files are used internally to verify that a given structure can be
regenerated from its internal coordinates.

Internal coordinates may also be exported as OpenSCAD data arrays for
generating 3D printed protein models.  OpenSCAD software is provided as
proof-of-concept for generating such models.

The following classes comprise the core functionality for processing internal
coordinates and are sufficiently related and coupled to place them together in
this module:

IC_Chain: Extends Biopython Chain on .internal_coord attribute.
    Manages connected sequence of residues and chain breaks; methods generally
    apply IC_Residue methods along chain.

IC_Residue: Extends for Biopython Residue on .internal_coord attribute.
    Most control and methods of interest are in this class, see API.

Dihedron: four joined atoms forming a dihedral angle.
    Dihedral angle, homogeneous atom coordinates in local coordinate space,
    references to relevant Hedra and IC_Residue.  Methods to compute
    residue dihedral angles, bond angles and bond lengths.

Hedron: three joined atoms forming a plane.
    Contains homogeneous atom coordinates in local coordinate space as well as
    bond lengths and angle between them.

Edron: base class for Hedron and Dihedron classes.
    Tuple of AtomKeys comprising child, string ID, mainchain membership boolean
    and other routines common for both Hedra and Dihedra.  Implements rich
    comparison.

AtomKey: keys (dictionary and string) for referencing atom sequences.
    Capture residue and disorder/occupancy information, provides a
    no-whitespace key for .pic files, and implements rich comparison.

Custom exception classes: HedronMatchError and MissingAtomError
"""
if TYPE_CHECKING:
    ...
HKT = Tuple["AtomKey", "AtomKey", "AtomKey"]
DKT = Tuple["AtomKey", "AtomKey", "AtomKey", "AtomKey"]
EKT = Union[HKT, DKT]
BKT = Tuple["AtomKey", "AtomKey"]
HACS = ...
DACS = Tuple[numpy.array, numpy.array, numpy.array, numpy.array]
class IC_Chain:
    """Class to extend Biopython Chain with internal coordinate data.

    Attributes
    ----------
    chain: biopython Chain object reference
        The Chain object this extends

    initNCaC: AtomKey indexed dictionary of N, Ca, C atom coordinates.
        NCaCKeys start chain segments (first residue or after chain break).
        These 3 atoms define the coordinate space for a contiguous chain segment,
        as initially specified by PDB or mmCIF file.

    MaxPeptideBond: **Class** attribute to detect chain breaks.
        Override for fully contiguous chains with some very long bonds - e.g.
        for 3D printing (OpenSCAD output) a structure with fully disordered
        (missing) residues.

    ordered_aa_ic_list: list of IC_Residue objects
        IC_Residue objects ic algorithms can process (e.g. no waters)

    hedra: dict indexed by 3-tuples of AtomKeys
        Hedra forming residues in this chain

    hedraLen: int length of hedra dict

    hedraNdx: dict mapping hedra AtomKeys to numpy array data

    dihedra: dict indexed by 4-tuples of AtomKeys
        Dihedra forming (overlapping) this residue

    dihedraLen: int length of dihedra dict

    dihedraNdx: dict mapping dihedra AtomKeys to numpy array data

    atomArray: numpy array of homogeneous atom coords for chain

    atomArrayIndex: dict mapping AtomKeys to atomArray indexes

    numpy arrays for vector processing of chain di/hedra:

    hedraIC: length-angle-length entries for each hedron

    hAtoms: homogeneous atom coordinates (3x4) of hedra, central atom at origin

    hAtomsR: hAtoms in reverse order

    hAtoms_needs_update: booleans indicating whether hAtoms represent hedraIC

    dihedraIC: dihedral angles for each dihedron

    dAtoms: homogeneous atom coordinates (4x4) of dihedra, second atom at origin

    dAtoms_needs_update: booleans indicating whether dAtoms represent dihedraIC

    Methods
    -------
    internal_to_atom_coordinates(verbose, start, fin)
        Process ic data to Residue/Atom coordinates; calls assemble_residues()
        followed by coords_to_structure()
    assemble_residues(verbose, start, fin)
        Generate IC_Residue atom coords from internal coordinates
    coords_to_structure()
        update Biopython Residue.Atom coords from IC_Residue coords for all
        Residues with IC_Residue attributes
    atom_to_internal_coordinates(verbose)
        Calculate dihedrals, angles, bond lengths (internal coordinates) for
        Atom data
    link_residues()
        Call link_dihedra() on each IC_Residue (needs rprev, rnext set)
    set_residues()
        Add .internal_coord attribute for all Residues in parent Chain, populate
        ordered_aa_ic_list, set IC_Residue rprev, rnext or initNCaC coordinates
    write_SCAD()
        Write OpenSCAD matrices for internal coordinate data comprising chain

    """
    MaxPeptideBond = ...
    def __init__(self, parent: Chain, verbose: bool = ...) -> None:
        """Initialize IC_Chain object, with or without residue/Atom data.

        :param parent: Biopython Chain object
            Chain object this extends
        """
        ...
    
    def clear_ic(self): # -> None:
        """Clear residue internal_coord settings for this chain."""
        ...
    
    def set_residues(self, verbose: bool = ...) -> None:
        """Initialize internal_coord data for loaded Residues.

        Add IC_Residue as .internal_coord attribute for each Residue in parent
        Chain; populate ordered_aa_ic_list with IC_Residue references for residues
        which can be built (amino acids and some hetatms); set rprev and rnext
        on each sequential IC_Residue, populate initNCaC at start and after
        chain breaks.
        """
        ...
    
    def link_residues(self) -> None:
        """link_dihedra() for each IC_Residue; needs rprev, rnext set.

        Called by PICIO:read_PIC() after finished reading chain
        """
        ...
    
    def assemble_residues(self, verbose: bool = ..., start: Optional[int] = ..., fin: Optional[int] = ...) -> None:
        """Generate IC_Residue atom coords from internal coordinates.

        Filter positions between start and fin if set, find appropriate start
        coordinates for each residue and pass to IC_Residue.assemble()

        :param verbose bool: default False
            describe runtime problems
        :param: start, fin lists
            sequence position, insert code for begin, end of subregion to
            process

        """
        ...
    
    def coords_to_structure(self) -> None:
        """Promote all ic atom_coords to Biopython Residue/Atom coords.

        IC atom_coords are homogeneous [4], Biopython atom coords are XYZ [3].
        """
        ...
    
    def init_edra(self) -> None:
        """Create chain level di/hedra arrays.

        If called by read_PIC, self.di/hedra = {} and object tree has IC data.
        -> build chain arrays from IC data

        If called at start of atom_to_internal_coords, self.di/hedra fully
        populated.  -> create empty chain numpy arrays

        In both cases, fix di/hedra object attributes to be views on
        chain-level array data
        """
        ...
    
    def init_atom_coords(self) -> None:
        """Set chain level di/hedra initial coord arrays from IC_Residue data."""
        ...
    
    def internal_to_atom_coordinates(self, verbose: bool = ..., start: Optional[int] = ..., fin: Optional[int] = ..., promote: Optional[bool] = ...) -> None:
        """Process, IC data to Residue/Atom coords.

        Not yet vectorized.

        :param verbose bool: default False
            describe runtime problems
        :param: start, fin lists
            sequence position, insert code for begin, end of subregion to
            process
        :param promote bool: default True
            If True (the default) copy result atom XYZ coordinates to
            Biopython Atom objects for access by other Biopython methods;
            otherwise, updated atom coordinates must be accessed through
            IC_Residue and hedron objects.
        """
        ...
    
    def atom_to_internal_coordinates(self, verbose: bool = ...) -> None:
        """Calculate dihedrals, angles, bond lengths for Atom data.

        :param verbose bool: default False
            describe runtime problems
        """
        ...
    
    def write_SCAD(self, fp: TextIO, backboneOnly: bool) -> None:
        """Write self to file fp as OpenSCAD data matrices.

        Works with write_SCAD() and embedded OpenSCAD routines in SCADIO.py.
        The OpenSCAD code explicitly creates spheres and cylinders to
        represent atoms and bonds in a 3D model.  Options are available
        to support rotatable bonds and magnetic hydrogen bonds.

        Matrices are written to link, enumerate and describe residues,
        dihedra, hedra, and chains, mirroring contents of the relevant IC_*
        data structures.

        The OpenSCAD matrix of hedra has additional information as follows:

        * the atom and bond state (single, double, resonance) are logged
          so that covalent radii may be used for atom spheres in the 3D models

        * bonds and atoms are tracked so that each is only created once

        * bond options for rotation and magnet holders for hydrogen bonds
          may be specified

        Note the application of IC_Chain attribute MaxPeptideBond: missing
        residues may be linked (joining chain segments with arbitrarily long
        bonds) by setting this to a large value.

        All ALTLOC (disordered) residues and atoms are written to the output model.
        """
        ...
    


class IC_Residue:
    """Class to extend Biopython Residue with internal coordinate data.

    Parameters
    ----------
    parent: biopython Residue object this class extends
    NO_ALTLOC: bool default False
    Disable processing of ALTLOC atoms if True, use only selected atoms.

    Attributes
    ----------
    residue: Biopython Residue object reference
        The Residue object this extends
    hedra: dict indexed by 3-tuples of AtomKeys
        Hedra forming this residue
    dihedra: dict indexed by 4-tuples of AtomKeys
        Dihedra forming (overlapping) this residue
    rprev, rnext: lists of IC_Residue objects
        References to adjacent (bonded, not missing, possibly disordered)
        residues in chain
    atom_coords: AtomKey indexed dict of numpy [4] arrays
        Local copy of atom homogeneous coordinates [4] for work
        distinct from Bopython Residue/Atom values
    alt_ids: list of char
        AltLoc IDs from PDB file
    bfactors: dict
        AtomKey indexed B-factors as read from PDB file
    NCaCKey: List of tuples of AtomKeys
        List of tuples of N, Ca, C backbone atom AtomKeys; usually only 1
        but more if backbone altlocs. Set by link_dihedra()
    is20AA: bool
        True if residue is one of 20 standard amino acids, based on
        Residue resname
    accept_atoms: tuple
        list of PDB atom names to use when generatiing internal coordinates.
        Default is:

        `accept_atoms = accept_mainchain + accept_hydrogens`

        to exclude hydrogens in internal coordinates and generated PDB files,
        override as:

        `IC_Residue.accept_atoms = IC_Residue.accept_mainchain`

        to get only backbone atoms plus amide proton, use:

        `IC_Residue.accept_atoms = IC_Residue.accept_mainchain + ('H',)`

        to convert D atoms to H, set `AtomKey.d2h = True` and use:

        `IC_Residue.accept_atoms = accept_mainchain + accept_hydrogens + accept_deuteriums`

        There is currently no option to output internal coordinates with D
        instead of H

    accept_resnames: tuple
        list of 3-letter residue names for HETATMs to accept when generating
        internal coordinates from atoms.  HETATM sidechain will be ignored, but normal
        backbone atoms (N, CA, C, O, CB) will be included.  Currently only
        CYG, YCM and UNK; override at your own risk.  To generate
        sidechain, add appropriate entries to ic_data_sidechains in
        ic_data.py and support in atom_to_internal_coordinates()
    gly_Cbeta: bool default False
        override class variable to True to generate internal coordinates for
        glycine CB atoms in atom_to_internal_coordinates().

        `IC_Residue.gly_Cbeta = True`
    allBonds: bool default False
        whereas a PDB file just specifies atoms, OpenSCAD output for 3D printing
        needs all bonds specified explicitly - otherwise, e.g. PHE rings will not
        be closed.  This variable is managed by the Write_SCAD() code and enables
        this.
    cic: IC_Chain default None
        parent chain IC_Chain object

    scale: optional float
        used for OpenSCAD output to generate gly_Cbeta bond length

    Methods
    -------
    applyMtx()
        multiply all IC_Residue atom_cords by passed matrix
    assemble(atomCoordsIn, resetLocation, verbose)
        Compute atom coordinates for this residue from internal coordinates
    atm241(coord)
        Convert 1x3 cartesian coords to 4x1 homogeneous coords
    coords_to_residue()
        Convert homogeneous atom_coords to Biopython cartesian Atom coords
    atom_to_internal_coordinates(verbose)
        Create hedra and dihedra for atom coordinates
    get_angle()
        Return angle for passed key
    get_length()
        Return bond length for specified pair
    link_dihedra()
        Link dihedra to this residue, form id3_dh_index
    load_PIC(edron)
        Process parsed (di-/h-)edron data from PIC file
    pick_angle()
        Find Hedron or Dihedron for passed key
    pick_length()
        Find hedra for passed AtomKey pair
    rak(atom info)
        Residue AtomKey - per residue AtomKey result cache
    set_angle()
        Set angle for passed key (no position updates)
    set_length()
        Set bond length in all relevant hedra for specified pair
    write_PIC(pdbid, chainId, s)
        Generate PIC format strings for this residue

    """
    accept_resnames = ...
    AllBonds: bool = ...
    def __init__(self, parent: Residue, NO_ALTLOC: bool = ...) -> None:
        """Initialize IC_Residue with parent Biopython Residue.

        :param parent: Biopython Residue object
            The Biopython Residue this object extends
        :param NO_ALTLOC: bool default False
            Option to disable processing altloc disordered atoms, use selected.
        """
        ...
    
    def rak(self, atm: Union[str, Atom]) -> AtomKey:
        """Cache calls to AtomKey for this residue."""
        ...
    
    def build_rak_cache(self) -> None:
        """Create explicit entries for for atoms so don't miss altlocs."""
        ...
    
    accept_mainchain = ...
    accept_hydrogens = ...
    accept_deuteriums = ...
    accept_atoms = ...
    gly_Cbeta = ...
    @staticmethod
    def atm241(coord: numpy.array) -> numpy.array:
        """Convert 1x3 cartesian coordinates to 4x1 homogeneous coordinates."""
        ...
    
    def __repr__(self) -> str:
        """Print string is parent Residue ID."""
        ...
    
    def pretty_str(self) -> str:
        """Nice string for residue ID."""
        ...
    
    def load_PIC(self, edron: Dict[str, str]) -> None:
        """Process parsed (di-/h-)edron data from PIC file.

        :param edron: parse dictionary from Edron.edron_re
        """
        ...
    
    def link_dihedra(self, verbose: bool = ...) -> None:
        """Housekeeping after loading all residues and dihedra.

        - Link dihedra to this residue
        - form id3_dh_index
        - form ak_set
        - set NCaCKey to be available AtomKeys
        """
        ...
    
    def set_flexible(self) -> None:
        """For OpenSCAD, mark N-CA and CA-C bonds to be flexible joints."""
        ...
    
    def set_hbond(self) -> None:
        """For OpenSCAD, mark H-N and C-O bonds to be hbonds (magnets)."""
        ...
    
    def default_startpos(self) -> Dict[AtomKey, numpy.array]:
        """Generate default N-Ca-C coordinates to build this residue from."""
        ...
    
    def get_startpos(self) -> Dict[AtomKey, numpy.array]:
        """Find N-Ca-C coordinates to build this residue from."""
        ...
    
    def clear_transforms(self): # -> None:
        """Set cst and rcst attributes to none before assemble()."""
        ...
    
    def assemble(self, resetLocation: bool = ..., verbose: bool = ...) -> Union[Dict[AtomKey, numpy.array], Dict[HKT, numpy.array], None]:
        """Compute atom coordinates for this residue from internal coordinates.

        Join dihedrons starting from N-CA-C and N-CA-CB hedrons, computing protein
        space coordinates for backbone and sidechain atoms

        Sets forward and reverse transforms on each Dihedron to convert from
        protein coordinates to dihedron space coordinates for first three
        atoms (see coord_space())

        Not vectorized (yet).

        **Algorithm**

        form double-ended queue, start with c-ca-n, o-c-ca, n-ca-cb, n-ca-c.

        if resetLocation=True, use initial coords from generating dihedron
        for n-ca-c initial positions (result in dihedron coordinate space)

        while queue not empty
            get 3-atom hedron key

            for each dihedron starting with hedron key (1st hedron of dihedron)

                if have coordinates for all 4 atoms already
                    add 2nd hedron key to back of queue
                else if have coordinates for 1st 3 atoms
                    compute forward and reverse transforms to take 1st 3 atoms
                    to/from dihedron initial coordinate space

                    use reverse transform to get position of 4th atom in
                    current coordinates from dihedron initial coordinates

                    add 2nd hedron key to back of queue
                else
                    ordering failed, put hedron key at back of queue and hope
                    next time we have 1st 3 atom positions (should not happen)

        loop terminates (queue drains) as hedron keys which do not start any
        dihedra are removed without action

        :param resetLocation: bool default False
            - Option to ignore start location and orient so N-Ca-C hedron
            at origin.

        :returns:
            Dict of AtomKey -> homogeneous atom coords for residue in protein space
            relative to previous residue

        """
        ...
    
    def atom_to_internal_coordinates(self, verbose: bool = ...) -> None:
        """Create hedra and dihedra for atom coordinates.

        :param verbose: bool default False
            warn about missing N, Ca, C backbone atoms.
        """
        ...
    
    def build_glyCB(self, gCBd: Dihedron): # -> None:
        """Populate values for Gly C-beta, rest of chain complete.

        Data averaged from Sep 2019 Dunbrack cullpdb_pc20_res2.2_R1.0
        restricted to structures with amide protons.

        Ala avg rotation of OCCACB from NCACO query:
        select avg(g.rslt) as avg_rslt, stddev(g.rslt) as sd_rslt, count(*)
        from
        (select f.d1d, f.d2d,
        (case when f.rslt > 0 then f.rslt-360.0 else f.rslt end) as rslt
        from (select d1.angle as d1d, d2.angle as d2d,
        (d2.angle - d1.angle) as rslt from dihedron d1,
        dihedron d2 where d1.rdh_class='AOACACAACB' and
        d2.rdh_class='ANACAACAO' and d1.pdb=d2.pdb and d1.chn=d2.chn
        and d1.res=d2.res) as f) as g

        | avg_rslt          | sd_rslt          | count   |
        | -122.682194862932 | 5.04403040513919 | 14098   |

        """
        ...
    
    def write_PIC(self, pdbid: str = ..., chainid: str = ..., s: str = ...) -> str:
        """Write PIC format lines for this residue.

        :param str pdbid: PDB idcode string; default 0PDB
        :param str chainid: PDB Chain ID character; default A
        :param str s: result string to add to
        """
        ...
    
    def coords_to_residue(self, rnext: bool = ...) -> None:
        """Convert self.atom_coords to biopython Residue Atom coords.

        Copy homogeneous IC_Residue atom_coords to self.residue cartesian
        Biopython Atom coords.

        :param rnext: bool default False
            next IC_Residue has no atom coords due to missing atoms, so try to
            populate with any available coords calculated for this residue
            di/hedra extending into next
        """
        ...
    
    relative_atom_re = ...
    def pick_angle(self, angle_key: Union[EKT, str]) -> Optional[Union[Hedron, Dihedron]]:
        """Get Hedron or Dihedron for angle_key.

        :param angle_key:
            - tuple of 3 or 4 AtomKeys
            - string of atom names ('CA') separated by :'s
            - string of [-1, 0, 1]<atom name> separated by ':'s. -1 is
              previous residue, 0 is this residue, 1 is next residue
            - psi, phi, omg, omega, chi1, chi2, chi3, chi4, chi5
            - tau (N-CA-C angle) see Richardson1981
            - except for tuples of AtomKeys, no option to access alternate disordered atoms

        Observe that a residue's phi and omega dihedrals, as well as the hedra
        comprising them (including the N:Ca:C tau hedron), are stored in the
        n-1 di/hedra sets; this is handled here, but may be an issue if accessing
        directly.

        The following are equivalent (except for sidechains with non-carbon
        atoms for chi2)::

            ric = r.internal_coord
            print(
                r,
                ric.get_angle("psi"),
                ric.get_angle("phi"),
                ric.get_angle("omg"),
                ric.get_angle("tau"),
                ric.get_angle("chi2"),
            )
            print(
                r,
                ric.get_angle("N:CA:C:1N"),
                ric.get_angle("-1C:N:CA:C"),
                ric.get_angle("-1CA:-1C:N:CA"),
                ric.get_angle("N:CA:C"),
                ric.get_angle("CA:CB:CG:CD"),
            )

        See ic_data.py for detail of atoms in the enumerated sidechain angles
        and the backbone angles which do not span the peptide bond. Using 's'
        for current residue ('self') and 'n' for next residue, the spanning
        angles are::

                (sN, sCA, sC, nN)   # psi
                (sCA, sC, nN, nCA)  # omega i+1
                (sC, nN, nCA, nC)   # phi i+1
                (sCA, sC, nN)
                (sC, nN, nCA)
                (nN, nCA, nC)       # tau i+1

        :return: Matching Hedron, Dihedron, or None.
        """
        ...
    
    def get_angle(self, angle_key: Union[EKT, str]) -> Optional[float]:
        """Get dihedron or hedron angle for specified key.

        See pick_angle() for key specifications.
        """
        ...
    
    def set_angle(self, angle_key: Union[EKT, str], v: float): # -> None:
        """Set dihedron or hedron angle for specified key.

        See pick_angle() for key specifications.
        """
        ...
    
    def pick_length(self, ak_spec: Union[str, BKT]) -> Tuple[Optional[List[Hedron]], Optional[BKT]]:
        """Get list of hedra containing specified atom pair.

        :param ak_spec:
            - tuple of two AtomKeys
            - string: two atom names separated by ':', e.g. 'N:CA' with
              optional position specifier relative to self, e.g. '-1C:N' for
              preceding peptide bond.

        The following are equivalent::

            ric = r.internal_coord
            print(
                r,
                ric.get_length("0C:1N"),
            )
            print(
                r,
                None
                if not ric.rnext
                else ric.get_length((ric.rak("C"), ric.rnext[0].rak("N"))),
            )

        :return: list of hedra containing specified atom pair, tuple of atom keys
        """
        ...
    
    def get_length(self, ak_spec: Union[str, BKT]) -> Optional[float]:
        """Get bond length for specified atom pair.

        See pick_length() for ak_spec.
        """
        ...
    
    def set_length(self, ak_spec: Union[str, BKT], val: float) -> None:
        """Set bond length for specified atom pair.

        See pick_length() for ak_spec.
        """
        ...
    
    def applyMtx(self, mtx: numpy.array) -> None:
        """Apply matrix to atom_coords for this IC_Residue."""
        ...
    


class Edron:
    """Base class for Hedron and Dihedron classes.

    Supports rich comparison based on lists of AtomKeys.

    Attributes
    ----------
    aks: tuple
        3 (hedron) or 4 (dihedron) AtomKeys defining this di/hedron
    id: str
        ':'-joined string of AtomKeys for this di/hedron
    needs_update: bool
        indicates di/hedron local atom_coords do NOT reflect current di/hedron
        angle and length values in hedron local coordinate space
    dh_class: str
        sequence of atoms (no position or residue) comprising di/hedron
        for statistics
    rdh_class: str
        sequence of residue, atoms comprising di/hedron for statistics
    edron_re: compiled regex (Class Attribute)
        A compiled regular expression matching string IDs for Hedron
        and Dihedron objects
    cic: IC_Chain reference
        Chain internal coords object containing this hedron

    Methods
    -------
    gen_key([AtomKey, ...] or AtomKey, ...) (Static Method)
        generate a ':'-joined string of AtomKey Ids
    gen_acs(atom_coords)
        generate tuple of atom coords for keys in self.aks
    is_backbone()
        Return True if all aks atoms are N, Ca, C or O

    """
    edron_re = ...
    @staticmethod
    def gen_key(lst: Union[List[str], List[AtomKey]]) -> str:
        """Generate string of ':'-joined AtomKey strings from input.

        :param lst: list of AtomKey objects or id strings
        """
        ...
    
    def __init__(self, *args: Union[List[AtomKey], EKT], **kwargs: str) -> None:
        """Initialize Edron with sequence of AtomKeys.

        Acceptable input:

            [ AtomKey, ... ]  : list of AtomKeys
            AtomKey, ...      : sequence of AtomKeys as args
            {'a1': str, 'a2': str, ... }  : dict of AtomKeys as 'a1', 'a2' ...
        """
        ...
    
    def gen_acs(self, atom_coords: Dict[AtomKey, numpy.array]) -> Tuple[numpy.array, ...]:
        """Generate tuple of atom coord arrays for keys in self.aks.

        :param atom_coords: AtomKey dict of atom coords for residue
        :raises: MissingAtomError any atoms in self.aks missing coordinates
        """
        ...
    
    def is_backbone(self) -> bool:
        """Report True for contains only N, C, CA, O, H atoms."""
        ...
    
    def __repr__(self) -> str:
        """Tuple of AtomKeys is default repr string."""
        ...
    
    def __hash__(self) -> int:
        """Hash calculated at init from aks tuple."""
        ...
    
    def __eq__(self, other: object) -> bool:
        """Test for equality."""
        ...
    
    def __ne__(self, other: object) -> bool:
        """Test for inequality."""
        ...
    
    def __gt__(self, other: object) -> bool:
        """Test greater than."""
        ...
    
    def __ge__(self, other: object) -> bool:
        """Test greater or equal."""
        ...
    
    def __lt__(self, other: object) -> bool:
        """Test less than."""
        ...
    
    def __le__(self, other: object) -> bool:
        """Test less or equal."""
        ...
    


class Hedron(Edron):
    """Class to represent three joined atoms forming a plane.

    Contains atom coordinates in local coordinate space, central atom
    at origin.  Stored in two orientations, with the 3rd (forward) or
    first (reversed) atom on the +Z axis.

    Attributes
    ----------
    lal: numpy array of len12, angle, len23
        len12 = distance between 1st and 2nd atom
        angle = angle (degrees) formed by 3 atoms
        len23 = distance between 2nd and 3rd atoms

    atoms: 3x4 numpy arrray (view on chain array)
        3 homogeneous atoms comprising hedron, 1st on XZ, 2nd at origin, 3rd on +Z
    atomsR: 3x4 numpy array (view on chain array)
        atoms reversed, 1st on +Z, 2nd at origin, 3rd on XZ plane

    Methods
    -------
    get_length()
        get bond length for specified atom pair
    set_length()
        set bond length for specified atom pair
    angle(), len12(), len23()
        gettters and setters for relevant attributes (angle in degrees)
    """
    def __init__(self, *args: Union[List[AtomKey], HKT], **kwargs: str) -> None:
        """Initialize Hedron with sequence of AtomKeys, kwargs.

        Acceptable input:
            As for Edron, plus optional 'len12', 'angle', 'len23'
            keyworded values.
        """
        ...
    
    def __repr__(self) -> str:
        """Print string for Hedron object."""
        ...
    
    @property
    def angle(self) -> float:
        """Get this hedron angle."""
        ...
    
    @angle.setter
    def angle(self, angle_deg) -> None:
        """Set this hedron angle; sets needs_update."""
        ...
    
    @property
    def len12(self): # -> Any | float:
        """Get first length for Hedron."""
        ...
    
    @len12.setter
    def len12(self, len): # -> None:
        """Set first length for Hedron; sets needs_update."""
        ...
    
    @property
    def len23(self) -> float:
        """Get second length for Hedron."""
        ...
    
    @len23.setter
    def len23(self, len): # -> None:
        """Set second length for Hedron; sets needs_update."""
        ...
    
    def get_length(self, ak_tpl: BKT) -> Optional[float]:
        """Get bond length for specified atom pair.

        :param ak_tpl: tuple of AtomKeys
            pair of atoms in this Hedron
        """
        ...
    
    def set_length(self, ak_tpl: BKT, newLength: float): # -> None:
        """Set bond length for specified atom pair; sets needs_update.

        :param ak_tpl: tuple of AtomKeys
            pair of atoms in this Hedron
        """
        ...
    


class Dihedron(Edron):
    """Class to represent four joined atoms forming a dihedral angle.

    Attributes
    ----------
    angle: float
        Measurement or specification of dihedral angle in degrees
    hedron1, hedron2: Hedron object references
        The two hedra which form the dihedral angle
    h1key, h2key: tuples of AtomKeys
        Hash keys for hedron1 and hedron2
    id3,id32: tuples of AtomKeys
        First 3 and second 3 atoms comprising dihedron; hxkey orders may differ
    initial_coords: tuple[4] of numpy arrays [4]
        Local atom coords for 4 atoms, [0] on XZ plane, [1] at origin,
        [2] on +Z, [3] rotated by dihedral
    a4_pre_rotation: numpy array [4]
        4th atom of dihedral aligned to XZ plane (angle not applied)
    ic_residue: IC_Residue object reference
        IC_Residue object containing this dihedral
    reverse: bool
        Indicates order of atoms in dihedron is reversed from order of atoms
        in hedra (configured by set_hedra())
    cst, rcst: numpy array [4][4]
        transforms to and from coordinate space defined by first hedron.
        set by IC_Residue.assemble().  defined by id3 order NOT h1key order
        (atoms may be reversed between these two)

    Methods
    -------
    set_hedra()
        work out hedra keys and orientation for this dihedron
    angle()
        getter/setter for dihdral angle in degrees

    """
    def __init__(self, *args: Union[List[AtomKey], DKT], **kwargs: str) -> None:
        """Initialize Dihedron with sequence of AtomKeys and optional dihedral angle.

        Acceptable input:
            As for Edron, plus optional 'dihedral' keyworded angle value.
        """
        ...
    
    def __repr__(self) -> str:
        """Print string for Dihedron object."""
        ...
    
    def set_hedra(self) -> Tuple[bool, Hedron, Hedron]:
        """Work out hedra keys and set rev flag."""
        ...
    
    @property
    def angle(self) -> float:
        """Get dihedral angle."""
        ...
    
    @angle.setter
    def angle(self, dangle_deg: float) -> None:
        """Save new dihedral angle; sets needs_update.

        faster to modify IC_Chain level arrays directly.

        N.B. dihedron (i-1)C-N-CA-CB is ignored if O exists.
        C-beta is by default placed using O-C-CA-CB, but O is missing
        in some PDB file residues, which means the sidechain cannot be
        placed.  The alternate CB path (i-1)C-N-CA-CB is provided to
        circumvent this, but if this is needed then it must be adjusted in
        conjunction with PHI ((i-1)C-N-CA-C) as they overlap.

        :param dangle_deg: float new dihedral angle in degrees
        """
        ...
    


class AtomKey:
    """Class for dict keys to reference atom coordinates.

    AtomKeys capture residue and disorder information together, and
    provide a no-whitespace string key for .pic files.

    Supports rich comparison and multiple ways to instantiate.

    AtomKeys contain:
     residue position, insertion code, 1 or 3 char residue name,
     atom name, altloc, and occupancy

    Attributes
    ----------
    akl: tuple
        All six fields of AtomKey
    fieldNames: tuple (Class Attribute)
        Mapping of key index positions to names
    fields: namedtuple (Class Attribute)
        Mapping of field names to index positions
    id: str
        '_'-joined AtomKey fields, excluding 'None' fields
    atom_re: compiled regex (Class Attribute)
        A compiled regular expression matching the string form of the key
    endnum_re: compiled regex (Class Attribute)
        A compiled regular expresion capturing digits at end of a string
    d2h: bool (Class Attribute)
        Convert D atoms to H on input; must also modify IC_Residue.accept_atoms
    missing: bool default False
        AtomKey __init__'d from string is probably missing, set this flag to
        note the issue (not set here)

    Methods
    -------
    altloc_match(other)
        Returns True if this AtomKey matches other AtomKey excluding altloc
        and occupancy fields

    """
    atom_re = ...
    endnum_re = ...
    fieldNames = ...
    fieldsDef = ...
    fields = ...
    d2h = ...
    def __init__(self, *args: Union[IC_Residue, Atom, List, Dict, str], **kwargs: str) -> None:
        """Initialize AtomKey with residue and atom data.

        Examples of acceptable input:
            (<IC_Residue>, 'CA', ...)    : IC_Residue with atom info
            (<IC_Residue>, <Atom>)       : IC_Residue with Biopython Atom
            ([52, None, 'G', 'CA', ...])  : list of ordered data fields
            (52, None, 'G', 'CA', ...)    : multiple ordered arguments
            ({respos: 52, icode: None, atm: 'CA', ...}) : dict with fieldNames
            (respos: 52, icode: None, atm: 'CA', ...) : kwargs with fieldNames
            52_G_CA, 52B_G_CA, 52_G_CA_0.33, 52_G_CA_B_0.33  : id strings
        """
        ...
    
    def __repr__(self) -> str:
        """Repr string from id."""
        ...
    
    def __hash__(self) -> int:
        """Hash calculated at init from akl tuple."""
        ...
    
    _backbone_sort_keys = ...
    _sidechain_sort_keys = ...
    _greek_sort_keys = ...
    def altloc_match(self, other: AtomKey) -> bool:
        """Test AtomKey match other discounting occupancy and altloc."""
        ...
    
    def __ne__(self, other: object) -> bool:
        """Test for inequality."""
        ...
    
    def __eq__(self, other: object) -> bool:
        """Test for equality."""
        ...
    
    def __gt__(self, other: object) -> bool:
        """Test greater than."""
        ...
    
    def __ge__(self, other: object) -> bool:
        """Test greater or equal."""
        ...
    
    def __lt__(self, other: object) -> bool:
        """Test less than."""
        ...
    
    def __le__(self, other: object) -> bool:
        """Test less or equal."""
        ...
    


def set_accuracy_95(num: float) -> float:
    """Reduce floating point accuracy to 9.5 (xxxx.xxxxx).

    Used by Hedron and Dihedron classes writing PIC and SCAD files.
    :param float num: input number
    :returns: float with specified accuracy
    """
    ...

class HedronMatchError(Exception):
    """Cannot find hedron in residue for given key."""
    ...


class MissingAtomError(Exception):
    """Missing atom coordinates for hedron or dihedron."""
    ...


