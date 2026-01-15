"""
Region classes for the cxregions package.

This module contains all region-related classes including the base JuliaRegion class
and specific region types like Interior1CRegion, Exterior1CRegion, and Annulus.
"""

import juliacall
import numpy as np

# Import the Julia module setup from the main package
from . import jl, JLCR


def wrap_jl_region(jul):
    """Wrap a Julia region object in the appropriate Python class.
    
    Parameters
    ----------
    jul : juliacall.AnyValue
        A Julia region object from ComplexRegions.jl
        
    Returns
    -------
    JuliaRegion
        The appropriate Python region wrapper
        
    Raises
    ------
    ValueError
        If the argument is not a Julia object or not a recognized region type
    """
    if not isinstance(jul, juliacall.AnyValue):  # type: ignore
        raise ValueError("Argument to wrap_jl_region is not a Julia object")
    if jl.isa(jul, JLCR.AbstractRegion):
        if jl.isa(jul, JLCR.Annulus):
            return Annulus(jul)
        elif jl.isa(jul, JLCR.ExteriorSimplyConnectedRegion):
            return Exterior1CRegion(jul)
        elif jl.isa(jul, JLCR.InteriorSimplyConnectedRegion):
            return Interior1CRegion(jul)
        elif jl.isa(jul, JLCR.ExteriorRegion):
            return ExteriorRegion(jul)
        elif jl.isa(jul, JLCR.InteriorConnectedRegion):
            return InteriorConnectedRegion(jul)
        else:
            raise ValueError(f"Julia type {jl.typeof(jul)} not recognized for wrapping")
    else:
        raise ValueError(f"Julia type {jl.typeof(jul)} not recognized for wrapping")


class JuliaRegion:
    """Base class for wrapping Julia region objects from ComplexRegions.jl.
    
    This class provides a Python interface to Julia region objects, which represent
    areas in the complex plane bounded by curves.
    
    Parameters
    ----------
    julia_obj : juliacall.AnyValue
        A Julia region object from ComplexRegions.jl
        
    Attributes
    ----------
    julia : juliacall.AnyValue
        The underlying Julia region object
    """
    
    def __init__(self, julia_obj):
        """Initialize a JuliaRegion wrapper.
        
        Parameters
        ----------
        julia_obj : juliacall.AnyValue
            A Julia region object from ComplexRegions.jl
            
        Raises
        ------
        ValueError
            If julia_obj is not a valid Julia region object
        """
        if isinstance(julia_obj, juliacall.AnyValue):  # type: ignore
            if jl.isa(julia_obj, JLCR.AbstractRegion):
                self.julia = julia_obj
        else:
            raise ValueError("Invalid argument to Region constructor")
        
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

    def contains(self, z=None):
        """Check if a point is contained in the region.
        
        Parameters
        ----------
        z : complex, optional
            Point to test for containment
            
        Returns
        -------
        bool
            True if z is in the region, False otherwise
        """
        if z is not None:
            return getattr(JLCR, "in")(z, self.julia)
        else:
            getattr(JLCR, "in")(self.julia)

    def innerboundary(self):
        """Get the inner boundary curves of the region.
        
        Returns
        -------
        JuliaPath or list of JuliaPath
            Inner boundary curve(s) of the region
        """
        from .paths import JuliaPath  # Import here to avoid circular imports
        b = JLCR.innerboundary(self.julia)
        if isinstance(b, juliacall.VectorValue):  # type: ignore
            return [JuliaPath(j) for j in b]
        else:
            return JuliaPath(b)

    def outerboundary(self):
        """Get the outer boundary curves of the region.
        
        Returns
        -------
        JuliaPath or list of JuliaPath
            Outer boundary curve(s) of the region
        """
        from .paths import JuliaPath  # Import here to avoid circular imports
        b = JLCR.outerboundary(self.julia)
        if isinstance(b, juliacall.VectorValue):  # type: ignore
            paths = []
            for j in b:
                paths.append(JuliaPath(j))
            return paths
        else:
            return JuliaPath(b)

    def union(self, other):
        """Compute the union of this region with another.
        
        Parameters
        ----------
        other : JuliaRegion
            Another region to union with
            
        Returns
        -------
        JuliaRegion
            Union of the two regions
        """
        r = JLCR.union(self.julia, other.julia)
        return JuliaRegion(r)
    
    def intersect(self, other):
        """Compute the intersection of this region with another.
        
        Parameters
        ----------
        other : JuliaRegion
            Another region to intersect with
            
        Returns
        -------
        JuliaRegion
            Intersection of the two regions
        """
        r = JLCR.intersect(self.julia, other.julia)
        return JuliaRegion(r)


def Jordan(c):
    """Construct a Jordan curve from a ClosedPath or ClosedCurve object.
    
    Parameters
    ----------
    c : ClosedPath, ClosedCurve, or juliacall.AnyValue
        A closed curve or path object
        
    Returns
    -------
    ClosedPath or ClosedCurve
        The input object if already a Jordan curve, or wrapped appropriately
        
    Raises
    ------
    ValueError
        If the argument is not a valid closed curve or path
    """
    from .paths import ClosedPath, wrap_jl_path  # Import here to avoid circular imports
    from .curves import ClosedCurve, wrap_jl_curve
    
    if isinstance(c, ClosedPath) or isinstance(c, ClosedCurve):
        return c
    elif isinstance(c, juliacall.AnyValue):  # type: ignore
        if jl.isa(c, JLCR.AbstractClosedPath):
            return wrap_jl_path(c)
        elif jl.isa(c, JLCR.AbstractClosedCurve):
            return wrap_jl_curve(c)
        else:
            raise ValueError("Julia argument to Jordan not recognized")
    else:
        raise ValueError("Argument to Jordan not recognized")


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
    from .paths import JuliaPath
    
    if isinstance(p, JuliaCurve) or isinstance(p, JuliaPath):
        return p.julia
    else:
        return p


class Exterior1CRegion(JuliaRegion):
    """Simply connected exterior region.
    
    This represents the region outside a single closed curve, extending to infinity.
    
    Parameters
    ----------
    boundary : ClosedCurve, ClosedPath, or juliacall.AnyValue
        The boundary curve of the region
        
    Attributes
    ----------
    boundary : ClosedCurve or ClosedPath
        The boundary curve of the region
        
    Examples
    --------
    >>> from cxregions.curves import Circle
    >>> circle = Circle(0, 1)
    >>> exterior = Exterior1CRegion(circle)
    """
    
    def __init__(self, boundary):
        if isinstance(boundary, juliacall.AnyValue):  # type: ignore
            if jl.isa(boundary, JLCR.ExteriorSimplyConnectedRegion) :
                self.julia = boundary
            else:
                raise ValueError("Invalid argument to Exterior1CRegion constructor")
        else:
            self.julia = JLCR.ExteriorSimplyConnectedRegion(boundary.julia)

        self.boundary = Jordan(JuliaRegion.get(self, "boundary"))

    def isfinite(self):
        """Check if the region is finite.
        
        Returns
        -------
        bool
            True if the boundary is finite, False otherwise
        """
        return self.boundary.isfinite()
    
    def __repr__(self):
        return str(f"Exterior simply connected region")


class ExteriorRegion(JuliaRegion):
    """Exterior region with multiple inner boundaries.
    
    This represents the region outside multiple closed curves, extending to infinity.
    
    Parameters
    ----------
    inner : list of ClosedCurve/ClosedPath or juliacall.AnyValue
        The inner boundary curves
        
    Attributes
    ----------
    inner : list of ClosedCurve/ClosedPath
        The inner boundary curves
        
    Examples
    --------
    >>> from cxregions.curves import Circle
    >>> circle1 = Circle(0, 1)
    >>> circle2 = Circle(3, 0.5)
    >>> exterior = ExteriorRegion([circle1, circle2])
    """
    
    def __init__(self, inner):
        if isinstance(inner, juliacall.AnyValue):  # type: ignore
            if jl.isa(inner, JLCR.ExteriorRegion):
                self.julia = inner
            else:
                raise ValueError("Invalid argument to ExteriorRegion constructor")
        else:
            innerb = juliacall.convert(jl.Vector, [get_julia(b) for b in inner])
            self.julia = JLCR.ExteriorRegion(innerb)

        b = JuliaRegion.get(self, "inner")
        self.inner = [Jordan(j) for j in b]

    def isfinite(self):
        """Check if the region is finite.
        
        Returns
        -------
        bool
            Always returns False for exterior regions
        """
        return False
    
    def __repr__(self):
        return f"Exterior region with {len(self.inner)} inner boundaries"


class Interior1CRegion(JuliaRegion):
    """Simply connected interior region.
    
    This represents the region inside a single closed curve.
    
    Parameters
    ----------
    boundary : ClosedCurve, ClosedPath, or juliacall.AnyValue
        The boundary curve of the region
        
    Attributes
    ----------
    boundary : ClosedCurve or ClosedPath
        The boundary curve of the region
        
    Examples
    --------
    >>> from cxregions.curves import Circle
    >>> circle = Circle(0, 1)
    >>> interior = Interior1CRegion(circle)
    """
    
    def __init__(self, boundary):
        if isinstance(boundary, juliacall.AnyValue):  # type: ignore
            if jl.isa(boundary, JLCR.InteriorSimplyConnectedRegion) :
                self.julia = boundary
            else:
                raise ValueError("Invalid argument to InteriorConnectedRegion constructor")
        else:
            self.julia = JLCR.InteriorSimplyConnectedRegion(boundary.julia)

        self.boundary = Jordan(JuliaRegion.get(self, "boundary"))

    def isfinite(self):
        """Check if the region is finite.
        
        Returns
        -------
        bool
            True if the boundary is finite, False otherwise
        """
        return self.boundary.isfinite()
    
    def __repr__(self):
        return str(f"Interior simply connected region")


class InteriorConnectedRegion(JuliaRegion):
    """Multiply connected interior region.
    
    This represents the region inside an outer boundary but outside inner boundaries.
    
    Parameters
    ----------
    outer : ClosedCurve, ClosedPath, or juliacall.AnyValue
        The outer boundary curve
    inner : list of ClosedCurve/ClosedPath, optional
        The inner boundary curves (holes)
        
    Attributes
    ----------
    outer : ClosedCurve or ClosedPath
        The outer boundary curve
    inner : list of ClosedCurve/ClosedPath
        The inner boundary curves
        
    Examples
    --------
    >>> from cxregions.curves import Circle
    >>> outer_circle = Circle(0, 2)
    >>> inner_circle = Circle(0, 1)
    >>> region = InteriorConnectedRegion(outer_circle, [inner_circle])
    """
    
    def __init__(self, outer, inner=[]):
        if isinstance(outer, juliacall.AnyValue):  # type: ignore
            if jl.isa(outer, JLCR.InteriorConnectedRegion) :
                self.julia = outer
            else:
                raise ValueError("Invalid argument to InteriorConnectedRegion constructor")
        else:
            innerb = juliacall.convert(jl.Vector, [get_julia(b) for b in inner])
            outerb = get_julia(outer)
            self.julia = JLCR.InteriorConnectedRegion(outerb, innerb)

        b = JuliaRegion.get(self, "inner")
        self.inner = [Jordan(j) for j in b]
        self.outer = JuliaRegion.get(self, "outer")

    def isfinite(self):
        """Check if the region is finite.
        
        Returns
        -------
        bool
            True if all boundaries are finite, False otherwise
        """
        return self.outer.isfinite() & all([b.isfinite() for b in self.inner])
    
    def __repr__(self):
        N = len(self.inner)
        return f"Interior {N+1}-connnected region"


class Annulus(InteriorConnectedRegion):
    """An annulus (ring-shaped region) between two circles.
    
    This represents the region between an inner and outer circle.
    
    Parameters
    ----------
    outer : float, Circle, or juliacall.AnyValue
        Outer radius, outer circle, or Julia Annulus object
    inner : float or Circle, optional
        Inner radius or inner circle
    center : complex, optional
        Center point (if constructing from radii), default is 0
        
    Attributes
    ----------
    inner : Circle
        Inner boundary circle
    outer : Circle
        Outer boundary circle
        
    Examples
    --------
    >>> # Annulus from radii
    >>> annulus1 = Annulus(2, 1, center=0)
    >>> # Annulus from circles
    >>> from cxregions.curves import Circle
    >>> inner_circle = Circle(0, 1)
    >>> outer_circle = Circle(0, 2)
    >>> annulus2 = Annulus(outer_circle, inner_circle)
    """
    
    def __init__(self, outer, inner=None, center=0j):
        from .curves import Circle  # Import here to avoid circular imports
        
        if isinstance(outer, juliacall.AnyValue):  # type: ignore
            if jl.isa(outer, JLCR.Annulus):
                self.julia = outer
            elif jl.isa(outer, JLCR.Circle) and jl.isa(inner, JLCR.Circle):
                self.julia = JLCR.Annulus(outer, inner)
            else:
                raise ValueError("Invalid argument to Annulus constructor")
        elif isinstance(inner, Circle) and isinstance(outer, Circle):
            self.julia = JLCR.Annulus(outer.julia, inner.julia)
        else:
            self.julia = JLCR.Annulus(outer, inner, center)

        self.inner = Circle(JuliaRegion.get(self, "inner"))
        self.outer = Circle(JuliaRegion.get(self, "outer"))

    def modulus(self):
        """Compute the modulus of the annulus.
        
        Returns
        -------
        float
            The modulus (related to the ratio of radii)
        """
        return JLCR.modulus(self.julia)
    
    def isfinite(self):
        """Check if the annulus is finite.
        
        Returns
        -------
        bool
            Always returns True for annuli
        """
        return True

    def __repr__(self):
        return f"Annulus centered at {self.inner.center} with radii {self.inner.radius} and {self.outer.radius}"


# Region construction functions
def between(curve1, curve2):
    """Construct the region between two closed curves.
    
    Parameters
    ----------
    curve1 : ClosedCurve or ClosedPath
        First boundary curve
    curve2 : ClosedCurve or ClosedPath
        Second boundary curve
        
    Returns
    -------
    InteriorConnectedRegion
        Region between the two curves
        
    Examples
    --------
    >>> from cxregions.curves import Circle
    >>> outer = Circle(0, 2)
    >>> inner = Circle(0, 1)
    >>> region = between(outer, inner)
    """
    r = JLCR.between(curve1.julia, curve2.julia)
    return InteriorConnectedRegion(r)


def interior(curve):
    """Construct the interior region of a closed curve.
    
    Parameters
    ----------
    curve : ClosedCurve or ClosedPath
        Boundary curve
        
    Returns
    -------
    Interior1CRegion
        Interior region of the curve
        
    Examples
    --------
    >>> from cxregions.curves import Circle
    >>> circle = Circle(0, 1)
    >>> region = interior(circle)
    """
    r = JLCR.interior(curve.julia)
    return Interior1CRegion(r)


def exterior(curve):
    """Construct the exterior region of a closed curve.
    
    Parameters
    ----------
    curve : ClosedCurve or ClosedPath
        Boundary curve
        
    Returns
    -------
    Exterior1CRegion
        Exterior region of the curve
        
    Examples
    --------
    >>> from cxregions.curves import Circle
    >>> circle = Circle(0, 1)
    >>> region = exterior(circle)
    """
    r = JLCR.exterior(curve.julia)
    return Exterior1CRegion(r)


def disk(center, radius):
    """Construct a disk as an interior region.
    
    Parameters
    ----------
    center : complex
        Center of the disk
    radius : float
        Radius of the disk
        
    Returns
    -------
    Interior1CRegion
        Disk region
        
    Examples
    --------
    >>> disk_region = disk(1+1j, 2)
    """
    r = JLCR.disk(center, radius)
    return Interior1CRegion(r)


def halfplane(l):
    """Construct a half-plane as an interior region from a Line.
    
    Parameters
    ----------
    l : Line
        Line that forms the boundary of the half-plane
        
    Returns
    -------
    Interior1CRegion
        Half-plane region
        
    Examples
    --------
    >>> from cxregions.curves import Line
    >>> line = Line(0, direction=1)
    >>> hp = halfplane(line)
    """
    r = JLCR.halfplane(l.julia)
    return Interior1CRegion(r)


# Pre-defined half-plane regions
from .curves import Line, Circle

upperhalfplane = halfplane(Line(0.0, direction=1.0))
"""Upper half-plane region (Im(z) > 0)."""

lowerhalfplane = halfplane(Line(0.0, direction=-1.0))
"""Lower half-plane region (Im(z) < 0)."""

lefthalfplane = halfplane(Line(0.0, direction=1.0j))
"""Left half-plane region (Re(z) < 0)."""

righthalfplane = halfplane(Line(0.0, direction=-1.0j))
"""Right half-plane region (Re(z) > 0)."""

unitdisk = interior(Circle(0.0, 1.0))
"""Unit disk region (|z| < 1)."""