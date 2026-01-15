"""
Path classes for the cxregions package.

This module contains all path-related classes including the base JuliaPath class
and specific path types like Polygon, CircularPolygon, and Rectangle.
"""

import juliacall
import numpy as np

# Import the Julia module setup from the main package
from . import jl, JLCR


def wrap_jl_path(jul):
    """Wrap a Julia path object in the appropriate Python class.
    
    Parameters
    ----------
    jul : juliacall.AnyValue
        A Julia path object from ComplexRegions.jl
        
    Returns
    -------
    Path
        The appropriate Python path wrapper
        
    Raises
    ------
    ValueError
        If the argument is not a Julia object or not a recognized path type
    """
    if not isinstance(jul, juliacall.AnyValue):  # type: ignore
        raise ValueError("Argument to wrap_jl_path is not a Julia object")
    if jl.isa(jul, JLCR.AbstractPath):
        if jl.isa(jul, JLCR.CircularPolygon):
            return CircularPolygon(jul)
        elif jl.isa(jul, JLCR.Polygon):
            return Polygon(jul)
        elif jl.isa(jul, JLCR.Rectangle):
            return Rectangle(jul)
        elif jl.isa(jul, JLCR.AbstractClosedPath):
            return ClosedPath(jul)
        else:
            return Path(jul)
    else:
        raise ValueError(f"Julia type {jl.typeof(jul)} not recognized for wrapping")


class JuliaPath:
    """Base class for wrapping Julia path objects from ComplexRegions.jl.
    
    This class provides a Python interface to Julia path objects, which are
    sequences of connected curves forming a continuous path.
    
    Parameters
    ----------
    julia_obj : juliacall.AnyValue
        A Julia path object from ComplexRegions.jl
        
    Attributes
    ----------
    julia : juliacall.AnyValue
        The underlying Julia path object
    """
    
    def __init__(self, julia_obj):
        """Initialize a JuliaPath wrapper.
        
        Parameters
        ----------
        julia_obj : juliacall.AnyValue
            A Julia path object from ComplexRegions.jl
            
        Raises
        ------
        ValueError
            If julia_obj is not a valid Julia path object
        """
        if isinstance(julia_obj, juliacall.AnyValue):  # type: ignore
            if jl.isa(julia_obj, JLCR.AbstractPath):
                self.julia = julia_obj
        else:
            raise ValueError("Invalid argument to Path constructor")

    def get(self, field):
        """Get a field from the underlying Julia object.
        
        Parameters
        ----------
        field : str
            Name of the field to retrieve
            
        Returns
        -------
        Any
            The value of the requested field
        """
        return jl.getproperty(self.julia, jl.Symbol(field))

    def length(self):
        """Get the number of curves in the path.
        
        Returns
        -------
        int
            Number of curves in the path
        """
        return JLCR.length(self.julia)
    
    def curves(self):
        """Get all curves that make up the path.
        
        Returns
        -------
        list of JuliaCurve
            List of curve objects that form the path
        """
        from .curves import Circle, Arc, Line, Segment, Ray, JuliaCurve
        
        curves = []
        for j in JLCR.curves(self.julia):
            if jl.isa(j, jl.Circle):
                curves.append(Circle(j))
            elif jl.isa(j, jl.Arc):
                curves.append(Arc(j))
            elif jl.isa(j, jl.Line):
                curves.append(Line(j))
            elif jl.isa(j, jl.Segment):
                curves.append(Segment(j))
            elif jl.isa(j, jl.Ray):
                curves.append(Ray(j))
            else:
                curves.append(JuliaCurve(j))
        return curves
    
    def curve(self, k):
        """Get the k-th curve in the path.
        
        Parameters
        ----------
        k : int
            Index of the curve (0-based)
            
        Returns
        -------
        JuliaCurve
            The k-th curve in the path
        """
        return self.curves()[k]

    def __getitem__(self, index):
        """Get a curve by index (enables path[i] syntax).
        
        Parameters
        ----------
        index : int
            Index of the curve
            
        Returns
        -------
        JuliaCurve
            The curve at the given index
        """
        return self.curve(index)

    def point(self, t):
        """Evaluate the path at parameter value t.
        
        Parameters
        ----------
        t : float
            Parameter value, typically in [0, 1]
            
        Returns
        -------
        complex
            Point on the path at parameter t
        """
        if isinstance(t, (list, tuple, np.ndarray)):
            # Use broadcasting in Julia for arrays of points
            res = jl.broadcast(self.julia, t)
            return np.array(res)
        else:
            return np.complex128(JLCR.point(self.julia, t))

    def __call__(self, t):
        """Evaluate the path at parameter value t (enables path(t) syntax).
        
        Parameters
        ----------
        t : float
            Parameter value, typically in [0, 1]
            
        Returns
        -------
        complex
            Point on the path at parameter t
        """
        return self.point(t)

    def arclength(self):
        """Compute the total arc length of the path.
        
        Returns
        -------
        float
            Total arc length of all curves in the path
        """
        return JLCR.arclength(self.julia)

    def tangent(self, t=0.):
        """Compute the tangent vector at parameter t.
        
        Parameters
        ----------
        t : float, optional
            Parameter value, default is 0
            
        Returns
        -------
        complex
            Tangent vector at parameter t
        """
        p = JLCR.tangent(self.julia, t)
        return np.complex128(p)

    def unittangent(self, t=0.):
        """Compute the unit tangent vector at parameter t.
        
        Parameters
        ----------
        t : float, optional
            Parameter value, default is 0
            
        Returns
        -------
        complex
            Unit tangent vector at parameter t
        """
        p = JLCR.unittangent(self.julia, t)
        return np.complex128(p)
    
    def normal(self, t=0.):
        """Compute the normal vector at parameter t.
        
        Parameters
        ----------
        t : float, optional
            Parameter value, default is 0
            
        Returns
        -------
        complex
            Normal vector at parameter t
        """
        p = JLCR.normal(self.julia, t)
        return np.complex128(p)

    def angles(self):
        """Get the turning angles at all vertices.
        
        Returns
        -------
        numpy.ndarray
            Array of turning angles at vertices
        """
        return np.array(JLCR.angles(self.julia))

    def vertices(self):
        """Get all vertices of the path.
        
        Returns
        -------
        numpy.ndarray
            Array of vertex coordinates
        """
        return np.array(JLCR.vertices(self.julia))
    
    def vertex(self, k):
        """Get the k-th vertex of the path.
        
        Parameters
        ----------
        k : int
            Index of the vertex
            
        Returns
        -------
        complex
            Coordinates of the k-th vertex
        """
        p = JLCR.vertex(self.julia, k)
        return np.complex128(p)

    def arg(self, z):
        """Find the parameter value corresponding to point z on the path.
        
        Parameters
        ----------
        z : complex
            Point to locate on the path
            
        Returns
        -------
        float or None
            Parameter value if z is on the path, None otherwise
        """
        return JLCR.arg(self.julia, z)

    def conj(self):
        """Return the complex conjugate of the path.
        
        Returns
        -------
        JuliaPath
            Complex conjugate of this path
        """
        p = JLCR.conj(self.julia)
        return type(self)(p)

    def reverse(self):
        """Return the path with reversed orientation.
        
        Returns
        -------
        JuliaPath
            Path with reversed orientation
        """
        p = JLCR.reverse(self.julia)
        return type(self)(p)
    
    def isfinite(self):
        """Check if the path has finite length.
        
        Returns
        -------
        bool
            True if the path is finite, False otherwise
        """
        return JLCR.isfinite(self.julia)
    
    def ispositive(self):
        """Check if the path has positive orientation.
        
        Returns
        -------
        bool
            True if positively oriented, False otherwise
        """
        return JLCR.ispositive(self.julia)

    def isreal(self):
        """Check if the path lies on the real axis.
        
        Returns
        -------
        bool
            True if the path is real, False otherwise
        """
        return JLCR.isreal(self.julia)

    def isapprox(self, other):
        """Check if this path is approximately equal to another.
        
        Parameters
        ----------
        other : JuliaPath
            Another path to compare with
            
        Returns
        -------
        bool
            True if paths are approximately equal, False otherwise
        """
        return JLCR.isapprox(self.julia, other.julia)

    def inv(self):
        """Compute the inversion of the path with respect to the unit circle.
        
        Returns
        -------
        JuliaPath
            Inverted path
        """
        p = JLCR.inv(self.julia)
        return type(self)(p)

    def reflect(self, z):
        """Reflect point z across the path.
        
        Parameters
        ----------
        z : complex
            Point to reflect
            
        Returns
        -------
        complex
            Reflected point
        """
        return JLCR.reflect(z, self.julia)

    def closest(self, z):
        """Find the closest point on the path to z.
        
        Parameters
        ----------
        z : complex
            Reference point
            
        Returns
        -------
        complex
            Closest point on the path to z
        """
        return JLCR.closest(z, self.julia)

    def dist(self, z):
        """Compute the distance from point z to the path.
        
        Parameters
        ----------
        z : complex
            Reference point
            
        Returns
        -------
        float
            Distance from z to the closest point on the path
        """
        return JLCR.dist(z, self.julia)
    
    def __add__(self, other):
        """Add a complex number to the path (translation).
        
        Parameters
        ----------
        other : complex
            Complex number to add
            
        Returns
        -------
        JuliaPath
            Translated path
        """
        julia_add = getattr(jl, "+")
        t = julia_add(self.julia, other)
        return type(self)(t)

    def __radd__(self, other):
        """Add a complex number to the path (translation).
        
        Parameters
        ----------
        other : complex
            Complex number to add
            
        Returns
        -------
        JuliaPath
            Translated path
        """
        julia_add = getattr(jl, "+")
        t = julia_add(other, self.julia)
        return type(self)(t)

    def __neg__(self):
        """Negate the path.
        
        Returns
        -------
        JuliaPath
            Negated path
        """
        julia_neg = getattr(jl, "-")
        t = julia_neg(self.julia)
        return type(self)(t)

    def __sub__(self, other):
        """Subtract a complex number from the path (translation).
        
        Parameters
        ----------
        other : complex
            Complex number to subtract
            
        Returns
        -------
        JuliaPath
            Translated path
        """
        julia_sub = getattr(jl, "-")
        t = julia_sub(self.julia, other)
        return type(self)(t)

    def __rsub__(self, other):
        """Subtract the path from a complex number.
        
        Parameters
        ----------
        other : complex
            Complex number
            
        Returns
        -------
        JuliaPath
            Resulting path
        """
        julia_sub = getattr(jl, "-")
        t = julia_sub(other, self.julia)
        return type(self)(t)
    
    def __mul__(self, other):
        """Multiply the path by a complex number (scaling and rotation).
        
        Parameters
        ----------
        other : complex
            Complex number to multiply by
            
        Returns
        -------
        JuliaPath
            Scaled and rotated path
        """
        julia_mul = getattr(jl, "*")
        t = julia_mul(self.julia, other)
        return type(self)(t)

    def __rmul__(self, other):
        """Multiply the path by a complex number (scaling and rotation).
        
        Parameters
        ----------
        other : complex
            Complex number to multiply by
            
        Returns
        -------
        JuliaPath
            Scaled and rotated path
        """
        julia_mul = getattr(jl, "*")
        t = julia_mul(other, self.julia)
        return type(self)(t)

    def __truediv__(self, other):
        """Divide the path by a complex number (scaling and rotation).
        
        Parameters
        ----------
        other : complex
            Complex number to divide by
            
        Returns
        -------
        JuliaPath
            Scaled and rotated path
        """
        julia_div = getattr(jl, "/")
        t = julia_div(self.julia, other)
        return type(self)(t)

    def intersect(self, other):
        """Find intersection points with another path.
        
        Parameters
        ----------
        other : JuliaPath
            Another path to intersect with
            
        Returns
        -------
        numpy.ndarray or other
            Array of intersection points or other result
        """
        z = JLCR.intersect(self.julia, other.julia)
        return z


class Path(JuliaPath):
    """A path composed of connected curves.
    
    This class represents a sequence of curves that form a continuous path.
    
    Parameters
    ----------
    curves : list of JuliaCurve or juliacall.AnyValue
        List of curves that form the path, or Julia Path object to wrap
        
    Examples
    --------
    >>> from cxregions.curves import Segment
    >>> seg1 = Segment(0, 1)
    >>> seg2 = Segment(1, 1+1j)
    >>> path = Path([seg1, seg2])
    """
    
    def __init__(self, curves):
        if isinstance(curves, juliacall.AnyValue):  # type: ignore
            if jl.isa(curves, JLCR.AbstractPath):
                self.julia = curves
            else:
                raise ValueError("Invalid argument to Path constructor")
        else:
            self.julia = JLCR.Path([c.julia for c in np.atleast_1d(curves)])
        
    def __repr__(self):
        N = len(self.curves())
        return f"Path with {N} curves"


class ClosedPath(Path):
    """A closed path that forms a Jordan curve.
    
    This class represents a closed sequence of curves that form a continuous
    closed path, enabling computation of winding numbers and containment tests.
    
    Parameters
    ----------
    curves : list of JuliaCurve, Path, or juliacall.AnyValue
        List of curves that form the closed path, or Julia ClosedPath object to wrap
        
    Examples
    --------
    >>> from cxregions.curves import Segment
    >>> seg1 = Segment(0, 1)
    >>> seg2 = Segment(1, 1+1j)
    >>> seg3 = Segment(1+1j, 1j)
    >>> seg4 = Segment(1j, 0)
    >>> closed_path = ClosedPath([seg1, seg2, seg3, seg4])
    """
    
    def __init__(self, curves):
        if isinstance(curves, juliacall.AnyValue):  # type: ignore
            if jl.isa(curves, JLCR.AbstractClosedPath):
                self.julia = curves
            else:
                raise ValueError("Invalid argument to ClosedPath constructor")
        elif isinstance(curves, Path):
            self.julia = JLCR.ClosedPath(curves.julia)
        else:
            self.julia = JLCR.ClosedPath([c.julia for c in np.atleast_1d(curves)])

    def winding(self, z):
        """Compute the winding number around point z.
        
        Parameters
        ----------
        z : complex
            Point around which to compute winding number
            
        Returns
        -------
        int
            Winding number (positive for counterclockwise orientation)
        """
        return JLCR.winding(self.julia, z)
            
    def isinside(self, z):
        """Check if point z is inside the closed path.
        
        Parameters
        ----------
        z : complex
            Point to test
            
        Returns
        -------
        bool
            True if z is inside the path, False otherwise
        """
        return JLCR.isinside(z, self.julia)

    def __repr__(self):
        N = len(self.curves())
        return f"Closed path with {N} curves"


def get_julia(p):
    """Extract the Julia object from a Python wrapper.
    
    Parameters
    ----------
    p : JuliaCurve, JuliaPath, or other
        Python wrapper object or raw Julia object
        
    Returns
    -------
    juliacall.AnyValue or other
        The underlying Julia object, or the input if not a wrapper
    """
    from .curves import JuliaCurve
    
    if isinstance(p, JuliaCurve) or isinstance(p, JuliaPath):
        return p.julia
    else:
        return p


class CircularPolygon(ClosedPath):
    """A polygon with circular arc sides.
    
    This class represents a polygon where each side can be either a straight
    line segment or a circular arc.
    
    Parameters
    ----------
    arg : list of curves or juliacall.AnyValue
        List of curves (segments and arcs) that form the polygon sides,
        or Julia CircularPolygon object to wrap
        
    Attributes
    ----------
    path : ClosedPath
        The underlying closed path
        
    Examples
    --------
    >>> from cxregions.curves import Arc, Segment
    >>> arc1 = Arc(1, 2 + 1j, 1j)
    >>> seg1 = Segment(1j, -1)
    >>> arc2 = Arc(-1, -0.5j, -1j)
    >>> seg2 = Segment(-1j, 1)
    >>> cpoly = CircularPolygon([arc1, seg1, arc2, seg2])
    """
    
    def __init__(self, arg):
        if isinstance(arg, juliacall.AnyValue):  # type: ignore
            if jl.isa(arg, JLCR.CircularPolygon):
                self.julia = arg
        else:
            vec = juliacall.convert(jl.Vector,[get_julia(a) for a in arg])
            self.julia = JLCR.CircularPolygon(vec)
        
        self.path = ClosedPath(JuliaPath.get(self, "path"))
    
    def sides(self):
        """Get all sides of the polygon.
        
        Returns
        -------
        list of JuliaCurve
            List of curves that form the polygon sides
        """
        return self.curves()
    
    def side(self, k):
        """Get the k-th side of the polygon.
        
        Parameters
        ----------
        k : int
            Index of the side
            
        Returns
        -------
        JuliaCurve
            The k-th side of the polygon
        """
        return self.curve(k)

    def __repr__(self):
        N = len(self.sides())
        return f"Circular polygon with {N} sides"


class Polygon(ClosedPath):
    """A polygon with straight sides.
    
    This class represents a polygon where all sides are straight line segments.
    
    Parameters
    ----------
    arg : list of complex or juliacall.AnyValue
        List of vertex coordinates or Julia Polygon object to wrap
        
    Attributes
    ----------
    path : ClosedPath
        The underlying closed path
        
    Examples
    --------
    >>> # Unit square
    >>> poly = Polygon([0, 1, 1+1j, 1j])
    >>> # More complex polygon
    >>> poly2 = Polygon([4, 4 + 3j, 3j, -2j, 6 - 2j, 6])
    """
    
    def __init__(self, arg):
        if isinstance(arg, juliacall.AnyValue):  # type: ignore
            if jl.isa(arg, JLCR.Polygon):
                self.julia = arg
        else:
            vec = juliacall.convert(jl.Vector,[get_julia(a) for a in arg])
            self.julia = JLCR.Polygon(vec)

        self.path = ClosedPath(JuliaPath.get(self, "path"))
    
    def sides(self):
        """Get all sides of the polygon.
        
        Returns
        -------
        list of JuliaCurve
            List of line segments that form the polygon sides
        """
        return self.curves()
    
    def side(self, k):
        """Get the k-th side of the polygon.
        
        Parameters
        ----------
        k : int
            Index of the side
            
        Returns
        -------
        JuliaCurve
            The k-th side of the polygon
        """
        return self.curve(k)
    
    def __repr__(self):
        N = len(self.sides())
        return f"Polygon with {N} sides"


class Rectangle(Polygon):
    """A rectangular polygon.
    
    This class represents a rectangle, which is a special case of a polygon
    with four sides at right angles.
    
    Parameters
    ----------
    a : complex, list, or juliacall.AnyValue
        Center point, list of vertices, or Julia Rectangle object to wrap
    b : array-like or complex, optional
        Radii array (if a is center) or opposite corner (if a is corner)
        
    Attributes
    ----------
    center : complex
        Center of the rectangle
    radii : array-like
        Half-widths in x and y directions
    rotation : float
        Rotation angle of the rectangle
    polygon : Polygon
        The underlying polygon representation
        
    Examples
    --------
    >>> # Rectangle from center and radii
    >>> rect1 = Rectangle(0+0j, np.array([1.0, 0.5]))
    >>> # Rectangle from opposite corners
    >>> rect2 = Rectangle(-1-1j, 1+1j)
    """
    
    def __init__(self, a, b=None):
        if isinstance(a, juliacall.AnyValue):  # type: ignore
            if jl.isa(a, JLCR.Rectangle):
                self.julia = a
            else:
                raise ValueError("Invalid argument to Rectangle constructor")
        else:
            if b is None:
                # hopefully, a vector of vertices was given
                self.julia = JLCR.rectangle(a)
            else:
                if np.ndim(a) == 0 and np.ndim(b) > 0:
                    # center and radii were given; use constructor
                    self.julia = JLCR.Rectangle(a, b)
                else:
                    # opposite corners were given; use rectangle function
                    self.julia = JLCR.rectangle(a, b)
        
        self.center = JuliaPath.get(self, "center")
        self.radii = JuliaPath.get(self, "radii")
        self.rotation = JuliaPath.get(self, "rotation")
        self.polygon = Polygon(JuliaPath.get(self, "polygon"))


# Utility functions
def n_gon(n):
    """Construct a regular n-gon as a Polygon object.
    
    Parameters
    ----------
    n : int
        Number of sides
        
    Returns
    -------
    Polygon
        Regular n-sided polygon centered at origin with unit circumradius
        
    Examples
    --------
    >>> # Regular hexagon
    >>> hex_poly = n_gon(6)
    >>> # Regular triangle
    >>> triangle = n_gon(3)
    """
    return Polygon(JLCR.n_gon(n))


def quad(rect):
    """Construct a quadrilateral region from a Rectangle.
    
    Parameters
    ----------
    rect : Rectangle
        Rectangle to convert to a region
        
    Returns
    -------
    Interior1CRegion
        Interior region of the rectangle
        
    Examples
    --------
    >>> rect = Rectangle(0+0j, np.array([1.0, 0.5]))
    >>> region = quad(rect)
    """
    from .regions import Interior1CRegion  # Import here to avoid circular imports
    r = JLCR.quad(rect.julia)
    return Interior1CRegion(r)