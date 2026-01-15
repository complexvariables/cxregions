"""
Curve classes for the cxregions package.

This module contains all curve-related classes including the base JuliaCurve class
and specific curve types like Line, Circle, Segment, Ray, and Arc.
"""

import juliacall
import numpy as np

# Import the Julia module setup from the main package
from . import jl, JLCR


def wrap_jl_curve(jul):
    """Wrap a Julia curve object in the appropriate Python class.
    
    Parameters
    ----------
    jul : juliacall.AnyValue
        A Julia curve object from ComplexRegions.jl
        
    Returns
    -------
    Curve
        The appropriate Python curve wrapper (Circle, Arc, Line, Segment, Ray, etc.)
        
    Raises
    ------
    ValueError
        If the argument is not a Julia object or not a recognized curve type
    """
    if not isinstance(jul, juliacall.AnyValue):  # type: ignore
        raise ValueError("Argument to wrap_jl_curve is not a Julia object")
    if jl.isa(jul, JLCR.AbstractCurve):
        if jl.isa(jul, JLCR.Circle):
            return Circle(jul)
        elif jl.isa(jul, JLCR.Arc):
            return Arc(jul)
        elif jl.isa(jul, JLCR.Line):
            return Line(jul)
        elif jl.isa(jul, JLCR.Segment):
            return Segment(jul)
        elif jl.isa(jul, JLCR.Ray):
            return Ray(jul)
        elif jl.isa(jul, JLCR.AbstractClosedCurve):
            return ClosedCurve(jul)
        else:
            return Curve(jul)
    else:
        raise ValueError(f"Julia type {jl.typeof(jul)} not recognized for wrapping")


class JuliaCurve:
    """Base class for wrapping Julia curve objects from ComplexRegions.jl.
    
    This class provides a Python interface to Julia curve objects, handling
    the conversion between Python and Julia types and providing access to
    geometric operations on curves.
    
    Parameters
    ----------
    julia_obj : juliacall.AnyValue
        A Julia curve object from ComplexRegions.jl
        
    Attributes
    ----------
    julia : juliacall.AnyValue
        The underlying Julia curve object
    """
    
    def __init__(self, julia_obj):
        """Initialize a JuliaCurve wrapper.
        
        Parameters
        ----------
        julia_obj : juliacall.AnyValue
            A Julia curve object from ComplexRegions.jl
            
        Raises
        ------
        ValueError
            If julia_obj is not a valid Julia curve object
        """
        if isinstance(julia_obj, juliacall.AnyValue):  # type: ignore
            if jl.isa(julia_obj, JLCR.AbstractCurve):
                self.julia = julia_obj
        else:
            raise ValueError("Invalid argument to Curve constructor")

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

    def point(self, t):
        """Evaluate the curve at parameter value t.
        
        Parameters
        ----------
        t : float
            Parameter value, typically in [0, 1]
            
        Returns
        -------
        complex
            Point on the curve at parameter t
        """
        if isinstance(t, (list, tuple, np.ndarray)):
            # Use broadcasting in Julia for arrays of points
            res = jl.broadcast(self.julia, t)
            return np.array(res)
        else:
            return np.complex128(JLCR.point(self.julia, t))

    def __call__(self, t):
        """Evaluate the curve at parameter value t (alias for point).
        
        Parameters
        ----------
        t : float
            Parameter value, typically in [0, 1]
        """
        return self.point(t)

    def arclength(self):
        """Compute the arc length of the curve.
        
        Returns
        -------
        float
            Total arc length of the curve
        """
        return JLCR.arclength(self.julia)

    def tangent(self, t=0.):
        """Compute the tangent at parameter t.
        
        Parameters
        ----------
        t : float, optional
            Parameter value, default is 0
            
        Returns
        -------
        complex
            Tangent at parameter t
        """
        p = JLCR.tangent(self.julia, t)
        return np.complex128(p)

    def unittangent(self, t=0.):
        """Compute the unit tangent at parameter t.
        
        Parameters
        ----------
        t : float, optional
            Parameter value, default is 0
            
        Returns
        -------
        complex
            Unit tangent at parameter t
        """
        p = JLCR.unittangent(self.julia, t)
        return np.complex128(p)
    
    def normal(self, t):
        """Compute the normal at parameter t.
        
        Parameters
        ----------
        t : float
            Parameter value
            
        Returns
        -------
        complex
            Normal at parameter t
        """
        p = JLCR.normal(self.julia, t)
        return np.complex128(p)

    def arg(self, z):
        """Find the parameter value corresponding to point z on the curve.
        
        Parameters
        ----------
        z : complex
            Point to locate on the curve
            
        Returns
        -------
        float or None
            Parameter value if z is on the curve, None otherwise
        """
        return JLCR.arg(self.julia, z)

    def conj(self):
        """Return the complex conjugate of the curve.
        
        Returns
        -------
        JuliaCurve
            Complex conjugate of this curve
        """
        c = JLCR.conj(self.julia)
        return type(self)(c)

    def reverse(self):
        """Return the curve with reversed orientation.
        
        Returns
        -------
        JuliaCurve
            Curve with reversed orientation
        """
        c = JLCR.reverse(self.julia)
        return type(self)(c)
    
    def isfinite(self):
        """Check if the curve has finite length.
        
        Returns
        -------
        bool
            True if the curve is finite, False otherwise
        """
        return JLCR.isfinite(self.julia)
    
    def ispositive(self):
        """Check if the curve has positive orientation.
        
        Returns
        -------
        bool
            True if positively oriented, False otherwise
        """
        return JLCR.ispositive(self.julia)

    def isreal(self):
        """Check if the curve lies on the real axis.
        
        Returns
        -------
        bool
            True if the curve is real, False otherwise
        """
        return JLCR.isreal(self.julia)

    def isapprox(self, other):
        """Check if this curve is approximately equal to another.
        
        Parameters
        ----------
        other : JuliaCurve
            Another curve to compare with
            
        Returns
        -------
        bool
            True if curves are approximately equal, False otherwise
        """
        return JLCR.isapprox(self.julia, other.julia)

    def inv(self):
        """Compute the inversion of the curve with respect to the origin.
        
        Returns
        -------
        juliacall.AnyValue
            Julia object representing the inverted curve
            
        Notes
        -----
        This method returns a raw Julia object. Subclasses should override
        this method to return properly wrapped Python objects.
        """
        # can't know the return type in general, so this must be wrapped by inheritors
        c = JLCR.inv(self.julia)
        return c

    def isleft(self, z):
        """Check if point z is to the left of the curve.
        
        Parameters
        ----------
        z : complex
            Point to test
            
        Returns
        -------
        bool
            True if z is to the left of the curve, False otherwise
        """
        return JLCR.isleft(z, self.julia)

    def isright(self, z):
        """Check if point z is to the right of the curve.
        
        Parameters
        ----------
        z : complex
            Point to test
            
        Returns
        -------
        bool
            True if z is to the right of the curve, False otherwise
        """
        return JLCR.isright(z, self.julia)
    
    def reflect(self, z):
        """Reflect point z across the curve.
        
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
        """Find the closest point on the curve to z.
        
        Parameters
        ----------
        z : complex
            Reference point
            
        Returns
        -------
        complex
            Closest point on the curve to z
        """
        return JLCR.closest(z, self.julia)

    def dist(self, z):
        """Compute the distance from point z to the curve.
        
        Parameters
        ----------
        z : complex
            Reference point
            
        Returns
        -------
        float
            Distance from z to the closest point on the curve
        """
        return JLCR.dist(z, self.julia)
    
    def __add__(self, other):
        """Add a complex number to the curve (translation).
        
        Parameters
        ----------
        other : complex
            Complex number to add
            
        Returns
        -------
        JuliaCurve
            Translated curve
        """
        julia_add = getattr(jl, "+")
        t = julia_add(self.julia, other)
        return type(self)(t)

    def __radd__(self, other):
        """Add a complex number to the curve (translation).
        
        Parameters
        ----------
        other : complex
            Complex number to add
            
        Returns
        -------
        JuliaCurve
            Translated curve
        """
        julia_add = getattr(jl, "+")
        t = julia_add(other, self.julia)
        return type(self)(t)

    def __neg__(self):
        """Negate the curve.
        
        Returns
        -------
        JuliaCurve
            Negated curve
        """
        julia_neg = getattr(jl, "-")
        t = julia_neg(self.julia)
        return type(self)(t)

    def __sub__(self, other):
        """Subtract a complex number from the curve (translation).
        
        Parameters
        ----------
        other : complex
            Complex number to subtract
            
        Returns
        -------
        JuliaCurve
            Translated curve
        """
        julia_sub = getattr(jl, "-")
        t = julia_sub(self.julia, other)
        return type(self)(t)

    def __rsub__(self, other):
        """Subtract the curve from a complex number.
        
        Parameters
        ----------
        other : complex
            Complex number
            
        Returns
        -------
        JuliaCurve
            Resulting curve
        """
        julia_sub = getattr(jl, "-")
        t = julia_sub(other, self.julia)
        return type(self)(t)
    
    def __mul__(self, other):
        """Multiply the curve by a complex number (scaling and rotation).
        
        Parameters
        ----------
        other : complex
            Complex number to multiply by
            
        Returns
        -------
        JuliaCurve
            Scaled and rotated curve
        """
        julia_mul = getattr(jl, "*")
        t = julia_mul(self.julia, other)
        return type(self)(t)

    def __rmul__(self, other):
        """Multiply the curve by a complex number (scaling and rotation).
        
        Parameters
        ----------
        other : complex
            Complex number to multiply by
            
        Returns
        -------
        JuliaCurve
            Scaled and rotated curve
        """
        julia_mul = getattr(jl, "*")
        t = julia_mul(other, self.julia)
        return type(self)(t)

    def __truediv__(self, other):
        """Divide the curve by a complex number (scaling and rotation).
        
        Parameters
        ----------
        other : complex
            Complex number to divide by
            
        Returns
        -------
        JuliaCurve
            Scaled and rotated curve
        """
        julia_div = getattr(jl, "/")
        t = julia_div(self.julia, other)
        return type(self)(t)

    def intersect(self, other):
        """Find intersection points with another curve.
        
        Parameters
        ----------
        other : JuliaCurve
            Another curve to intersect with
            
        Returns
        -------
        numpy.ndarray or JuliaCurve
            Array of intersection points, or a curve if the intersection
            is a continuous curve segment
        """
        z = JLCR.intersect(self.julia, other.julia)
        if isinstance(z, juliacall.VectorValue):  # type: ignore
            if len(z) == 0:
                return np.array([])
            else:
                return np.array(z)
        elif jl.isa(z, JLCR.AbstractCurve):
            return wrap_jl_curve(z)
        else:
            return z


class Curve(JuliaCurve):
    """A general parametric curve in the complex plane.
    
    This class represents a curve defined by a point function and optional
    tangent function over a parameter domain.
    
    Parameters
    ----------
    point : callable or juliacall.AnyValue
        Either a function that maps parameter values to complex points,
        or a Julia curve object to wrap
    tangent : callable
        Function that maps parameter values to tangent vectors
    domain : tuple of float, optional
        Parameter domain as (start, end), default is (0.0, 1.0)
        
    Examples
    --------
    >>> # Create a curve from a function
    >>> curve = Curve(lambda t: t + 1j*t**2, lambda t: 1 + 2j*t)
    """
    
    def __init__(self, point, tangent=None, domain=(0.0, 1.0)):
        if isinstance(point, juliacall.AnyValue):  # type: ignore
            if jl.isa(point, JLCR.Curve):
                self.julia = point
            else:
                raise ValueError("Invalid argument to Curve constructor")
        else:
            self.julia = JLCR.Curve(point, tangent, domain[0], domain[1])

    def inv(self):
        """Compute the inversion of the curve with respect to the origin.
        
        Returns
        -------
        Curve
            Inverted curve
        """
        c = JuliaCurve.inv(self)
        return type(self)(c)

    def __repr__(self):
        return str("Curve")


class ClosedCurve(Curve):
    """A closed parametric curve in the complex plane.
    
    This class represents a closed curve that forms a Jordan curve,
    enabling computation of winding numbers and interior/exterior tests.
    
    Parameters
    ----------
    point : callable or juliacall.AnyValue
        Either a function that maps parameter values to complex points,
        or a Julia closed curve object to wrap
    tangent : callable
        Function that maps parameter values to tangent vectors
    domain : tuple of float, optional
        Parameter domain as (start, end), default is (0.0, 1.0)
        
    Examples
    --------
    >>> # Create a closed curve (unit circle)
    >>> point = lambda t: np.exp(2j * np.pi * t)
    >>> tangent = lambda t: 2j * np.pi * np.exp(2j * np.pi * t)
    >>> curve = ClosedCurve(point, tangent)
    """
    
    def __init__(self, point, tangent=None, domain=(0.0, 1.0)):
        if isinstance(point, juliacall.AnyValue):  # type: ignore
            if jl.isa(point, JLCR.AbstractClosedCurve):
                self.julia = point
            else:
                raise ValueError("Invalid argument to ClosedCurve constructor")
        else:
            self.julia = JLCR.ClosedCurve(point, tangent, domain[0], domain[1])

    def winding(self, z):
        """Compute the winding number of the curve around point z.
        
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
    
    def hasinside(self, z):
        """Check if point z is in the interior of the closed curve.
        
        Parameters
        ----------
        z : complex
            Point to test

        Returns
        -------
        bool
            True if z is in the interior, False otherwise
        """
        return JLCR.isinside(z, self.julia)

    def __repr__(self):
        return str("Closed curve")


class Line(Curve):
    """An infinite straight line in the complex plane.
    
    A line can be constructed from two points, or from a point and a direction.
    Lines have infinite arc length and are always considered to have positive
    orientation.
    
    Parameters
    ----------
    a : complex or juliacall.AnyValue
        Either a point on the line, or a Julia Line object to wrap
    b : complex, optional
        Second point on the line (if constructing from two points)
    direction : complex, optional
        Direction vector of the line (if constructing from point and direction)
        
    Attributes
    ----------
    base : complex
        A base point on the line
    direction : complex
        Direction vector of the line
        
    Examples
    --------
    >>> # Line through two points
    >>> line1 = Line(0, 1+1j)
    >>> # Line through origin with given direction
    >>> line2 = Line(0, direction=1+1j)
    """
    
    def __init__(self, a, b=None, direction=None):
        if isinstance(a, juliacall.AnyValue): # type: ignore
            if jl.isa(a, JLCR.Line):
                self.julia = a
            else:
                raise ValueError("Invalid argument to Line constructor")
        elif b is not None:
            self.julia = JLCR.Line(a, b)
        else:
            self.julia = JLCR.Line(a, direction=direction)
        self.base = JuliaCurve.get(self, "base")
        self.direction = JuliaCurve.get(self, "direction")

    def arclength(self):
        """Return infinite arc length for lines.
        
        Returns
        -------
        float
            Always returns numpy.inf
        """
        return np.inf

    def ispositive(self):
        """Check if the line has positive orientation.
        
        Returns
        -------
        bool
            Always returns True for lines
        """
        return True

    def isfinite(self):
        """Check if the line is finite.
        
        Returns
        -------
        bool
            Always returns False for lines
        """
        return False

    def inv(self):
        """Compute the inversion of the line with respect to the origin.
        
        Returns
        -------
        Circle or Line
            Circle if the line doesn't pass through the origin, Line otherwise
        """
        c = JuliaCurve.inv(self)
        if jl.isa(c, JLCR.Circle):
            return Circle(c)
        else:
            return Line(c)

    def slope(self):
        """Get the slope of the line.
        
        Returns
        -------
        float
            Slope of the line (dy/dx)
        """
        return JLCR.slope(self.julia)

    def angle(self):
        """Get the angle of the line's direction vector.
        
        Returns
        -------
        float
            Angle in radians
        """
        return JLCR.angle(self.julia)
    
    def __repr__(self):
        return f"Line through {self.point(0.5)} at angle {self.angle() / np.pi} * pi"


class Circle(ClosedCurve):
    """A circle in the complex plane.
    
    Circles can be constructed in several ways:
    - From center and radius
    - From three points on the circle
    - From a Julia Circle object
    
    Parameters
    ----------
    a : complex or juliacall.AnyValue
        Center point, first point on circle, or Julia Circle object
    b : float or complex, optional
        Radius (if a is center) or second point on circle
    c : complex, optional
        Third point on circle (if constructing from three points)
    ccw : bool, optional
        Whether the circle has counterclockwise orientation, default True
        
    Attributes
    ----------
    center : complex
        Center of the circle
    radius : float
        Radius of the circle
    ccw : bool
        Whether the circle has counterclockwise orientation
        
    Examples
    --------
    >>> # Unit circle at origin
    >>> circle1 = Circle(0, 1)
    >>> # Circle from center and radius with clockwise orientation
    >>> circle2 = Circle(1+1j, 2, ccw=False)
    >>> # Circle through three points
    >>> circle3 = Circle(1, 1j, -1)
    """
    
    def __init__(self, a, b=None, c=None, ccw=True):
        if b is None:
            if isinstance(a, juliacall.AnyValue): # type: ignore
                if jl.isa(a, JLCR.Circle):
                    self.julia = a
            else:
                raise ValueError("Invalid argument to Circle constructor")
        elif c is None:
            self.julia = JLCR.Circle(a, b, ccw)
        else:
            self.julia = JLCR.Circle(a, b, c)
        
        self.radius = JuliaCurve.get(self, "radius")
        self.center = JuliaCurve.get(self, "center")
        self.ccw = JuliaCurve.get(self, "ccw")

    def ispositive(self):
        """Check if the circle has positive (counterclockwise) orientation.
        
        Returns
        -------
        bool
            True if counterclockwise, False if clockwise
        """
        return self.ccw

    def isfinite(self):
        """Check if the circle is finite.
        
        Returns
        -------
        bool
            Always returns True for circles
        """
        return True

    def inv(self):
        """Compute the inversion of the circle with respect to the origin.
        
        Returns
        -------
        Circle or Line
            Circle if the original circle doesn't pass through the origin,
            Line if it does
        """
        c = JuliaCurve.inv(self)
        if jl.isa(c, JLCR.Circle):
            return Circle(c)
        else:
            return Line(c)

    def __repr__(self):
        return f"Circle centered at {self.center} with radius {self.radius}"


class Segment(Curve):
    """A line segment in the complex plane.
    
    A segment is a finite portion of a line connecting two endpoints.
    
    Parameters
    ----------
    a : complex or juliacall.AnyValue
        First endpoint or Julia Segment object to wrap
    b : complex, optional
        Second endpoint (if a is first endpoint)
        
    Attributes
    ----------
    first : complex
        First endpoint of the segment
    last : complex
        Last endpoint of the segment
        
    Examples
    --------
    >>> # Segment from origin to (1,1)
    >>> seg = Segment(0, 1+1j)
    """
    
    def __init__(self, a, b=None):
        if isinstance(a, juliacall.AnyValue):  # type: ignore
            if jl.isa(a, JLCR.Segment):
                self.julia = a
            else:
                raise ValueError("Invalid argument to Segment constructor")
        else:
            self.julia = JLCR.Segment(a, b)
        
        self.first = JuliaCurve.get(self, "za")
        self.last = JuliaCurve.get(self, "zb")

    def inv(self):
        """Compute the inversion of the segment with respect to the origin.
        
        Returns
        -------
        Arc, Ray, or Segment
            The type depends on the geometry after inversion
        """
        c = JuliaCurve.inv(self)
        if jl.isa(c, JLCR.Arc):
            return Arc(c)
        elif jl.isa(c, JLCR.Ray):
            return Ray(c)
        else:
            return Segment(c)

    def __repr__(self):
        return f"Segment from {self.first} to {self.last}"


class Ray(Curve):
    """A semi-infinite ray in the complex plane.
    
    A ray starts at a base point and extends infinitely in a given direction.
    
    Parameters
    ----------
    base : complex or juliacall.AnyValue
        Starting point of the ray or Julia Ray object to wrap
    angle : float, optional
        Angle of the ray direction in radians
        
    Attributes
    ----------
    base : complex
        Starting point of the ray
    angle : float
        Angle of the ray direction
        
    Examples
    --------
    >>> # Ray from origin at 45 degrees
    >>> ray = Ray(0, np.pi/4)
    """
    
    def __init__(self, base, angle=None):
        if isinstance(base, juliacall.AnyValue):  # type: ignore
            if jl.isa(base, JLCR.Ray):
                self.julia = base
            else:
                raise ValueError("Invalid argument to Ray constructor")
        else:
            self.julia = JLCR.Ray(base, angle)
        self.base = JuliaCurve.get(self, "base")
        self.angle = JuliaCurve.get(self, "angle")

    def __repr__(self):
        return f"Ray from {self.base} at angle {self.angle / np.pi} * pi"


class Arc(Curve):
    """A circular arc in the complex plane.
    
    An arc is a portion of a circle between two points.
    
    Parameters
    ----------
    a : complex, Circle, or juliacall.AnyValue
        Start point, circle, or Julia Arc object to wrap
    b : complex, optional
        End point (if a is start point)
    c : complex, optional
        Center point (if constructing from start, end, center)
        
    Attributes
    ----------
    circle : Circle or Segment
        The underlying circle (or degenerate segment)
    start : float
        Starting parameter on the circle
    delta : float
        Parameter range of the arc
        
    Examples
    --------
    >>> # Quarter circle arc from 1 to i around origin
    >>> arc = Arc(1, 1j, 0)
    """
    
    def __init__(self, a, b=None, c=None):
        if isinstance(a, juliacall.AnyValue):  # type: ignore
            if jl.isa(a, JLCR.Arc):
                self.julia = a
            elif jl.isa(a, JLCR.Circle):
                self.julia = JLCR.Arc(a, b, c)
            else:
                raise ValueError("Invalid argument to Arc constructor")
        else:
            self.julia = JLCR.Arc(a, b, c)
        
        circ = JuliaCurve.get(self, "circle")
        try:
            self.circle = Circle(circ)
        except Exception:
            self.circle = Segment(circ)
        self.start = JuliaCurve.get(self, "start")
        self.delta = JuliaCurve.get(self, "delta")

    def inv(self):
        """Compute the inversion of the arc with respect to the origin.
        
        Returns
        -------
        Arc, Ray, or Segment
            The type depends on the geometry after inversion
        """
        c = JuliaCurve.inv(self)
        if jl.isa(c, JLCR.Arc):
            return Arc(c)
        elif jl.isa(c, JLCR.Ray):
            return Ray(c)
        else:
            return Segment(c)

    def __repr__(self):
        return f"Arc: fraction {self.delta} of {self.circle} from {self.start}"


# Module-level constants
unitcircle = Circle(0, 1)
"""Unit circle centered at the origin."""