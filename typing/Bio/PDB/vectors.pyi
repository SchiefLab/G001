"""
This type stub file was generated by pyright.
"""

import numpy
from typing import Optional, Tuple

"""Vector class, including rotation-related functions."""
def m2rotaxis(m): # -> tuple[float, Vector] | tuple[float | Any | Literal[0], Vector]:
    """Return angles, axis pair that corresponds to rotation matrix m.

    The case where ``m`` is the identity matrix corresponds to a singularity
    where any rotation axis is valid. In that case, ``Vector([1, 0, 0])``,
    is returned.
    """
    ...

def vector_to_axis(line, point):
    """Vector to axis method.

    Return the vector between a point and
    the closest point on a line (ie. the perpendicular
    projection of the point on the line).

    :type line: L{Vector}
    :param line: vector defining a line

    :type point: L{Vector}
    :param point: vector defining the point
    """
    ...

def rotaxis2m(theta, vector): # -> ndarray:
    """Calculate left multiplying rotation matrix.

    Calculate a left multiplying rotation matrix that rotates
    theta rad around vector.

    :type theta: float
    :param theta: the rotation angle

    :type vector: L{Vector}
    :param vector: the rotation axis

    :return: The rotation matrix, a 3x3 Numeric array.

    Examples
    --------
    >>> from numpy import pi
    >>> from Bio.PDB.vectors import rotaxis2m
    >>> from Bio.PDB.vectors import Vector
    >>> m = rotaxis2m(pi, Vector(1, 0, 0))
    >>> Vector(1, 2, 3).left_multiply(m)
    <Vector 1.00, -2.00, -3.00>

    """
    ...

rotaxis = ...
def refmat(p, q): # -> ndarray | Any:
    """Return a (left multiplying) matrix that mirrors p onto q.

    :type p,q: L{Vector}
    :return: The mirror operation, a 3x3 Numeric array.

    Examples
    --------
    >>> from Bio.PDB.vectors import refmat
    >>> p, q = Vector(1, 2, 3), Vector(2, 3, 5)
    >>> mirror = refmat(p, q)
    >>> qq = p.left_multiply(mirror)
    >>> print(q)
    <Vector 2.00, 3.00, 5.00>
    >>> print(qq)
    <Vector 1.21, 1.82, 3.03>

    """
    ...

def rotmat(p, q): # -> Any:
    """Return a (left multiplying) matrix that rotates p onto q.

    :param p: moving vector
    :type p: L{Vector}

    :param q: fixed vector
    :type q: L{Vector}

    :return: rotation matrix that rotates p onto q
    :rtype: 3x3 Numeric array

    Examples
    --------
    >>> from Bio.PDB.vectors import rotmat
    >>> p, q = Vector(1, 2, 3), Vector(2, 3, 5)
    >>> r = rotmat(p, q)
    >>> print(q)
    <Vector 2.00, 3.00, 5.00>
    >>> print(p)
    <Vector 1.00, 2.00, 3.00>
    >>> p.left_multiply(r)
    <Vector 1.21, 1.82, 3.03>

    """
    ...

def calc_angle(v1, v2, v3):
    """Calculate angle method.

    Calculate the angle between 3 vectors
    representing 3 connected points.

    :param v1, v2, v3: the tree points that define the angle
    :type v1, v2, v3: L{Vector}

    :return: angle
    :rtype: float
    """
    ...

def calc_dihedral(v1, v2, v3, v4):
    """Calculate dihedral angle method.

    Calculate the dihedral angle between 4 vectors
    representing 4 connected points. The angle is in
    ]-pi, pi].

    :param v1, v2, v3, v4: the four points that define the dihedral angle
    :type v1, v2, v3, v4: L{Vector}
    """
    ...

class Vector:
    """3D vector."""
    def __init__(self, x, y=..., z=...) -> None:
        """Initialize the class."""
        ...
    
    def __repr__(self): # -> str:
        """Return vector 3D coordinates."""
        ...
    
    def __neg__(self): # -> Vector:
        """Return Vector(-x, -y, -z)."""
        ...
    
    def __add__(self, other): # -> Vector:
        """Return Vector+other Vector or scalar."""
        ...
    
    def __sub__(self, other): # -> Vector:
        """Return Vector-other Vector or scalar."""
        ...
    
    def __mul__(self, other): # -> int:
        """Return Vector.Vector (dot product)."""
        ...
    
    def __truediv__(self, x): # -> Vector:
        """Return Vector(coords/a)."""
        ...
    
    def __pow__(self, other): # -> Vector:
        """Return VectorxVector (cross product) or Vectorxscalar."""
        ...
    
    def __getitem__(self, i): # -> Any:
        """Return value of array index i."""
        ...
    
    def __setitem__(self, i, value): # -> None:
        """Assign values to array index i."""
        ...
    
    def __contains__(self, i): # -> bool:
        """Validate if i is in array."""
        ...
    
    def norm(self): # -> Any:
        """Return vector norm."""
        ...
    
    def normsq(self): # -> int:
        """Return square of vector norm."""
        ...
    
    def normalize(self): # -> None:
        """Normalize the Vector object.

        Changes the state of ``self`` and doesn't return a value.
        If you need to chain function calls or create a new object
        use the ``normalized`` method.
        """
        ...
    
    def normalized(self): # -> Vector:
        """Return a normalized copy of the Vector.

        To avoid allocating new objects use the ``normalize`` method.
        """
        ...
    
    def angle(self, other): # -> Any:
        """Return angle between two vectors."""
        ...
    
    def get_array(self): # -> ndarray:
        """Return (a copy of) the array of coordinates."""
        ...
    
    def left_multiply(self, matrix): # -> Vector:
        """Return Vector=Matrix x Vector."""
        ...
    
    def right_multiply(self, matrix): # -> Vector:
        """Return Vector=Vector x Matrix."""
        ...
    
    def copy(self): # -> Vector:
        """Return a deep copy of the Vector."""
        ...
    


def homog_rot_mtx(angle_rads: float, axis: str) -> numpy.array:
    """Generate a 4x4 single-axis numpy rotation matrix.

    :param float angle_rads: the desired rotation angle in radians
    :param char axis: character specifying the rotation axis
    """
    ...

def set_Z_homog_rot_mtx(angle_rads: float, mtx: numpy.ndarray): # -> None:
    """Update existing Z rotation matrix to new angle."""
    ...

def set_Y_homog_rot_mtx(angle_rads: float, mtx: numpy.ndarray): # -> None:
    """Update existing Y rotation matrix to new angle."""
    ...

def set_X_homog_rot_mtx(angle_rads: float, mtx: numpy.ndarray): # -> None:
    """Update existing X rotation matrix to new angle."""
    ...

def homog_trans_mtx(x: float, y: float, z: float) -> numpy.array:
    """Generate a 4x4 numpy translation matrix.

    :param x, y, z: translation in each axis
    """
    ...

def set_homog_trans_mtx(x: float, y: float, z: float, mtx: numpy.ndarray): # -> None:
    """Update existing translation matrix to new values."""
    ...

def homog_scale_mtx(scale: float) -> numpy.array:
    """Generate a 4x4 numpy scaling matrix.

    :param float scale: scale multiplier
    """
    ...

def get_spherical_coordinates(xyz: numpy.array) -> Tuple[float, float, float]:
    """Compute spherical coordinates (r, azimuth, polar_angle) for X,Y,Z point.

    :param array xyz: column vector (3 row x 1 column numpy array)
    :return: tuple of r, azimuth, polar_angle for input coordinate
    """
    ...

gtm = ...
gmrz = ...
gmry = ...
gmrz2 = ...
def coord_space(a0: numpy.ndarray, a1: numpy.ndarray, a2: numpy.ndarray, rev: bool = ...) -> Tuple[numpy.ndarray, Optional[numpy.ndarray]]:
    """Generate transformation matrix to coordinate space defined by 3 points.

    New coordinate space will have:
        acs[0] on XZ plane
        acs[1] origin
        acs[2] on +Z axis

    :param numpy column array x3 acs: X,Y,Z column input coordinates x3
    :param bool rev: if True, also return reverse transformation matrix
        (to return from coord_space)
    :returns: 4x4 numpy array, x2 if rev=True
    """
    ...

def multi_rot_Z(angle_rads: numpy.ndarray) -> numpy.ndarray:
    """Create [entries] numpy Z rotation matrices for [entries] angles.

    :param entries: int number of matrices generated.
    :param angle_rads: numpy array of angles
    :returns: entries x 4 x 4 homogeneous rotation matrices
    """
    ...

def multi_rot_Y(angle_rads: numpy.ndarray) -> numpy.ndarray:
    """Create [entries] numpy Y rotation matrices for [entries] angles.

    :param entries: int number of matrices generated.
    :param angle_rads: numpy array of angles
    :returns: entries x 4 x 4 homogeneous rotation matrices
    """
    ...

