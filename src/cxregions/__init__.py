import juliacall
import numpy as np
jl = juliacall.newmodule("PyCR")
jl.seval('import Pkg')
installed = False
for v in jl.Pkg.dependencies().values():
    if v.name == "ComplexRegions":
        installed = True
        break
if not installed:
    jl.seval('Pkg.add("ComplexRegions")')
    
jl.seval("using ComplexRegions, PythonCall")

__all__ = ["Curve", "ClosedCurve", "Line", "Segment", "Circle", "Ray", "Arc", 
           "Path", "ClosedPath", "CircularPolygon", "Polygon", "Rectangle", "n_gon", "unitcircle",
           "ExteriorRegion", "InteriorRegion"]

class JuliaCurve:
    def __init__(self, julia_obj):
        if isinstance(julia_obj, juliacall.AnyValue):  # type: ignore
            if jl.isa(julia_obj, jl.ComplexRegions.AbstractCurve):
                self.julia = julia_obj
        else:
            raise ValueError("Invalid argument to Curve constructor")

    def get(self, field):
        return jl.getproperty(self.julia, jl.Symbol(field))

    def point(self, t):
        p = jl.ComplexRegions.point(self.julia, t)
        return np.complex128(p)

    def arclength(self):
        return jl.ComplexRegions.arclength(self.julia)

    def tangent(self, t=0.):
        p = jl.ComplexRegions.tangent(self.julia, t)
        return np.complex128(p)

    def unittangent(self, t=0.):
        p = jl.ComplexRegions.unittangent(self.julia, t)
        return np.complex128(p)
    
    def normal(self, t):
        p = jl.ComplexRegions.normal(self.julia, t)
        return np.complex128(p)

    def arg(self, z):
        return jl.ComplexRegions.arg(self.julia, z)

    def conj(self):
        c = jl.ComplexRegions.conj(self.julia)
        return type(self)(c)

    def reverse(self):
        c = jl.ComplexRegions.reverse(self.julia)
        return type(self)(c)
    
    def isfinite(self):
        return jl.ComplexRegions.isfinite(self.julia)
    
    def ispositive(self):
        return jl.ComplexRegions.ispositive(self.julia)

    def isreal(self):
        return jl.ComplexRegions.isreal(self.julia)

    def isapprox(self, other):
        return jl.ComplexRegions.isapprox(self.julia, other.julia)

    def inv(self):
        # can't know the return type in general, so this must be wrapped by inheritors
        c = jl.ComplexRegions.inv(self.julia)
        return c

    def isleft(self, z):
        return jl.ComplexRegions.isleft(z, self.julia)

    def isright(self, z):
        return jl.ComplexRegions.isright(z, self.julia)
    
    def reflect(self, z):
        return jl.ComplexRegions.reflect(z, self.julia)

    def closest(self, z):
        return jl.ComplexRegions.closest(z, self.julia)

    def dist(self, z):
        return jl.ComplexRegions.dist(z, self.julia)
    
    def __add__(self, other):
        julia_add = getattr(jl, "+")
        t = julia_add(self.julia, other)
        return type(self)(t)

    def __radd__(self, other):
        julia_add = getattr(jl, "+")
        t = julia_add(other, self.julia)
        return type(self)(t)

    def __neg__(self):
        julia_neg = getattr(jl, "-")
        t = julia_neg(self.julia)
        return type(self)(t)

    def __sub__(self, other):
        julia_sub = getattr(jl, "-")
        t = julia_sub(self.julia, other)
        return type(self)(t)

    def __rsub__(self, other):
        julia_sub = getattr(jl, "-")
        t = julia_sub(other, self.julia)
        return type(self)(t)
    
    def __mul__(self, other):
        julia_mul = getattr(jl, "*")
        t = julia_mul(self.julia, other)
        return type(self)(t)

    def __rmul__(self, other):
        julia_mul = getattr(jl, "*")
        t = julia_mul(other, self.julia)
        return type(self)(t)

    def __truediv__(self, other):
        julia_div = getattr(jl, "/")
        t = julia_div(self.julia, other)
        return type(self)(t)

    def intersect(self, other):
        z = jl.ComplexRegions.intersect(self.julia, other.julia)
        if isinstance(z, juliacall.VectorValue):  # type: ignore
            return np.array(z)
        elif jl.isa(z.julia, jl.Circle):
            return Circle(z.julia)
        elif jl.isa(z.julia, jl.Arc):
            return Arc(z.julia)
        elif jl.isa(z.julia, jl.Line):
            return Line(z.julia)
        elif jl.isa(z.julia, jl.Segment):
            return Segment(z.julia)
        elif jl.isa(z.julia, jl.Ray):
            return Ray(z.julia)
        else:
            return z

class Curve(JuliaCurve):
    def __init__(self, point, tangent=None, domain=(0.0, 1.0)):
        if isinstance(point, juliacall.AnyValue):  # type: ignore
            if jl.isa(point, jl.ComplexRegions.Curve):
                self.julia = point
            else:
                raise ValueError("Invalid argument to Curve constructor")
        else:
            self.julia = jl.ComplexRegions.Curve(point, tangent, domain[0], domain[1])

    def inv(self):
        c = JuliaCurve.inv(self)
        return type(self)(c)

    def __repr__(self):
        return str("Curve")

class ClosedCurve(Curve):
    def __init__(self, point, tangent=None, domain=(0.0, 1.0)):
        if isinstance(point, juliacall.AnyValue):  # type: ignore
            if jl.isa(point, jl.ComplexRegions.ClosedCurve):
                self.julia = point
            else:
                raise ValueError("Invalid argument to ClosedCurve constructor")
        else:
            self.julia = jl.ComplexRegions.ClosedCurve(point, tangent, domain[0], domain[1])

    def winding(self, z):
        return jl.ComplexRegions.winding(self.julia, z)

    def __repr__(self):
        return str("Closed curve")

class Line(Curve):
    def __init__(self, a, b=None, direction=None):
        if isinstance(a, juliacall.AnyValue): # type: ignore
            if jl.isa(a, jl.ComplexRegions.Line):
                self.julia = a
            else:
                raise ValueError("Invalid argument to Line constructor")
        elif b is not None:
            self.julia = jl.ComplexRegions.Line(a, b)
        else:
            self.julia = jl.ComplexRegions.Line(a, direction=direction)
        self.base = JuliaCurve.get(self, "base")
        self.direction = JuliaCurve.get(self, "direction")

    def arclength(self):
        return np.inf

    def ispositive(self):
        return True

    def isfinite(self):
        return False

    def inv(self):
        c = JuliaCurve.inv(self)
        if jl.isa(c, jl.ComplexRegions.Circle):
            return Circle(c)
        else:
            return Line(c)

    def slope(self):
        return jl.ComplexRegions.slope(self.julia)

    def angle(self):
        return jl.ComplexRegions.angle(self.julia)
    
    def __repr__(self):
        return f"Line through {self.point(0.5)} at angle {self.angle() / np.pi} * pi"

class Circle(ClosedCurve):
    def __init__(self, a, b=None, c=None, ccw=True):
        if b is None:
            if isinstance(a, juliacall.AnyValue): # type: ignore
                if jl.isa(a, jl.ComplexRegions.Circle):
                    self.julia = a
            else:
                raise ValueError("Invalid argument to Circle constructor")
        elif c is None:
            self.julia = jl.ComplexRegions.Circle(a, b, ccw)
        else:
            self.julia = jl.ComplexRegions.Circle(a, b, c)
        
        self.radius = JuliaCurve.get(self, "radius")
        self.center = JuliaCurve.get(self, "center")
        self.ccw = JuliaCurve.get(self, "ccw")

    def ispositive(self):
        return self.ccw()

    def isfinite(self):
        return True

    def inv(self):
        c = JuliaCurve.inv(self)
        if jl.isa(c, jl.ComplexRegions.Circle):
            return Circle(c)
        else:
            return Line(c)

    def __repr__(self):
        return f"Circle centered at {self.center} with radius {self.radius}"
 
class Segment(Curve):
    def __init__(self, a, b=None):
        if isinstance(a, juliacall.AnyValue):  # type: ignore
            if jl.isa(a, jl.ComplexRegions.Segment):
                self.julia = a
            else:
                raise ValueError("Invalid argument to Segment constructor")
        else:
            self.julia = jl.ComplexRegions.Segment(a, b)
        
        self.first = JuliaCurve.get(self, "za")
        self.last = JuliaCurve.get(self, "zb")

    def inv(self):
        c = JuliaCurve.inv(self)
        if jl.isa(c, jl.ComplexRegions.Arc):
            return Arc(c)
        elif jl.isa(c, jl.ComplexRegions.Ray):
            return Ray(c)
        else:
            return Segment(c)

    def __repr__(self):
        return f"Segment from {self.first} to {self.last}"

class Ray(Curve):
    def __init__(self, base, angle=None):
        if isinstance(base, juliacall.AnyValue):  # type: ignore
            if jl.isa(base, jl.ComplexRegions.Ray):
                self.julia = base
            else:
                raise ValueError("Invalid argument to Ray constructor")
        else:
            self.julia = jl.ComplexRegions.Ray(base, angle)
        self.base = JuliaCurve.get(self, "base")
        self.angle = JuliaCurve.get(self, "angle")

    def __repr__(self):
        return f"Ray from {self.base} at angle {self.angle / np.pi} * pi"

class Arc(Curve):
    def __init__(self, a, b=None, c=None):
        if isinstance(a, juliacall.AnyValue):  # type: ignore
            if jl.isa(a, jl.ComplexRegions.Arc):
                self.julia = a
            elif jl.isa(a, jl.ComplexRegions.Circle):
                self.julia = jl.ComplexRegions.Arc(a, b, c)
            else:
                raise ValueError("Invalid argument to Arc constructor")
        else:
            self.julia = jl.ComplexRegions.Arc(a, b, c)
        
        circ = JuliaCurve.get(self, "circle")
        try:
            self.circle = Circle(circ)
        except Exception:
            self.circle = Segment(circ)
        self.start = JuliaCurve.get(self, "start")
        self.delta = JuliaCurve.get(self, "delta")

    def inv(self):
        c = JuliaCurve.inv(self)
        if jl.isa(c, jl.ComplexRegions.Arc):
            return Arc(c)
        elif jl.isa(c, jl.ComplexRegions.Ray):
            return Ray(c)
        else:
            return Segment(c)

    def __repr__(self):
        return f"Arc: fraction {self.delta} of {self.circle} from {self.start}"
    
class JuliaPath:
    def __init__(self, julia_obj):
        if isinstance(julia_obj, juliacall.AnyValue):  # type: ignore
            if jl.isa(julia_obj, jl.ComplexRegions.AbstractPath):
                self.julia = julia_obj
        else:
            raise ValueError("Invalid argument to Path constructor")

    def get(self, field):
        return jl.getproperty(self.julia, jl.Symbol(field))

    def length(self):
        return jl.ComplexRegions.length(self.julia)
    
    def curves(self):
        curves = []
        for j in jl.ComplexRegions.curves(self.julia):
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
        return self.curves()[k]

    def __getitem__(self, index):
        return self.curve(index)

    def point(self, t):
        p = jl.ComplexRegions.point(self.julia, t)
        return np.complex128(p)

    def arclength(self):
        return jl.ComplexRegions.arclength(self.julia)

    def tangent(self, t=0.):
        p = jl.ComplexRegions.tangent(self.julia, t)
        return np.complex128(p)

    def unittangent(self, t=0.):
        p = jl.ComplexRegions.unittangent(self.julia, t)
        return np.complex128(p)
    
    def normal(self, t=0.):
        p = jl.ComplexRegions.normal(self.julia, t)
        return np.complex128(p)

    def angles(self):
        return np.array(jl.ComplexRegions.angles(self.julia))

    def vertices(self):
        return np.array(jl.ComplexRegions.vertices(self.julia))
    
    def vertex(self, k):
        p = jl.ComplexRegions.vertex(self.julia, k)
        return np.complex128(p)

    def arg(self, z):
        return jl.ComplexRegions.arg(self.julia, z)

    def conj(self):
        p = jl.ComplexRegions.conj(self.julia)
        return type(self)(p)

    def reverse(self):
        p = jl.ComplexRegions.reverse(self.julia)
        return type(self)(p)
    
    def isfinite(self):
        return jl.ComplexRegions.isfinite(self.julia)
    
    def ispositive(self):
        return jl.ComplexRegions.ispositive(self.julia)

    def isreal(self):
        return jl.ComplexRegions.isreal(self.julia)

    def isapprox(self, other):
        return jl.ComplexRegions.isapprox(self.julia, other.julia)

    def inv(self):
        p = jl.ComplexRegions.inv(self.julia)
        return type(self)(p)

    def reflect(self, z):
        return jl.ComplexRegions.reflect(z, self.julia)

    def closest(self, z):
        return jl.ComplexRegions.closest(z, self.julia)

    def dist(self, z):
        return jl.ComplexRegions.dist(z, self.julia)
    
    def __add__(self, other):
        julia_add = getattr(jl, "+")
        t = julia_add(self.julia, other)
        return type(self)(t)

    def __radd__(self, other):
        julia_add = getattr(jl, "+")
        t = julia_add(other, self.julia)
        return type(self)(t)

    def __neg__(self):
        julia_neg = getattr(jl, "-")
        t = julia_neg(self.julia)
        return type(self)(t)

    def __sub__(self, other):
        julia_sub = getattr(jl, "-")
        t = julia_sub(self.julia, other)
        return type(self)(t)

    def __rsub__(self, other):
        julia_sub = getattr(jl, "-")
        t = julia_sub(other, self.julia)
        return type(self)(t)
    
    def __mul__(self, other):
        julia_mul = getattr(jl, "*")
        t = julia_mul(self.julia, other)
        return type(self)(t)

    def __rmul__(self, other):
        julia_mul = getattr(jl, "*")
        t = julia_mul(other, self.julia)
        return type(self)(t)

    def __truediv__(self, other):
        julia_div = getattr(jl, "/")
        t = julia_div(self.julia, other)
        return type(self)(t)

    def intersect(self, other):
        z = jl.ComplexRegions.intersect(self.julia, other.julia)
        return z

class Path(JuliaPath):
    def __init__(self, curves):
        if isinstance(curves, juliacall.AnyValue):  # type: ignore
            if jl.isa(curves, jl.ComplexRegions.Path):
                self.julia = curves
            else:
                raise ValueError("Invalid argument to Path constructor")
        else:
            self.julia = jl.ComplexRegions.Path([c.julia for c in curves])
        
        self.curve = self.get("curve")
        
    def __repr__(self):
        N = len(self.curves())
        return f"Path with {N} curves"

class ClosedPath(Path):
    def __init__(self, curves):
        if isinstance(curves, juliacall.AnyValue):  # type: ignore
            if jl.isa(curves, jl.ComplexRegions.ClosedPath):
                self.julia = curves
            else:
                raise ValueError("Invalid argument to ClosedPath constructor")
        elif isinstance(curves, Path):
            self.julia = jl.ComplexRegions.ClosedPath(curves.julia)
        else:
            self.julia = jl.ComplexRegions.ClosedPath([c.julia for c in curves])
        
        self.curve = self.get("curve")

    def winding(self, z):
        return jl.ComplexRegions.winding(self.julia, z)
            
    def isinside(self, z):
        return jl.ComplexRegions.isinside(z, self.julia)

    def __repr__(self):
        N = len(self.curves())
        return f"Closed path with {N} curves"

def get_julia(p):
    if isinstance(p, JuliaCurve) or isinstance(p, JuliaPath):
        return p.julia
    else:
        return p

class CircularPolygon(ClosedPath):
    def __init__(self, arg):
        if isinstance(arg, juliacall.AnyValue):  # type: ignore
            if jl.isa(arg, jl.ComplexRegions.CircularPolygon):
                self.julia = arg
        else:
            vec = juliacall.convert(jl.Vector,[get_julia(a) for a in arg])
            self.julia = jl.ComplexRegions.CircularPolygon(vec)
        
        self.path = ClosedPath(JuliaPath.get(self, "path"))
    
    def sides(self):
        return self.curves()
    
    def side(self, k):
        return self.curve(k)

class Polygon(ClosedPath):
    def __init__(self, arg):
        if isinstance(arg, juliacall.AnyValue):  # type: ignore
            if jl.isa(arg, jl.ComplexRegions.Polygon):
                self.julia = arg
        else:
            vec = juliacall.convert(jl.Vector,[get_julia(a) for a in arg])
            self.julia = jl.ComplexRegions.Polygon(vec)

        self.path = ClosedPath(JuliaPath.get(self, "path"))
    
    def sides(self):
        return self.curves()
    
    def side(self, k):
        return self.curve(k)
    
class Rectangle(Polygon):
    def __init__(self, a, b=None):
        if isinstance(a, juliacall.AnyValue):  # type: ignore
            if jl.isa(a, jl.ComplexRegions.Rectangle):
                self.julia = a
            else:
                raise ValueError("Invalid argument to Rectangle constructor")
        else:
            if b is None:
                self.julia = jl.ComplexRegions.Rectangle(a)
            else:
                self.julia = jl.ComplexRegions.Rectangle(a, b)
        
        self.center = JuliaPath.get(self, "center")
        self.radii = JuliaPath.get(self, "radii")
        self.rotation = JuliaPath.get(self, "rotation")
        self.polygon = Polygon(JuliaPath.get(self, "polygon"))

##
unitcircle = Circle(0, 1)
def n_gon(n):
    """Construct a regular n-gon as a Polygon object."""
    return Polygon(jl.ComplexRegions.n_gon(n))

class JuliaRegion:
    def __init__(self, julia_obj):
        if isinstance(julia_obj, juliacall.AnyValue):  # type: ignore
            if jl.isa(julia_obj, jl.ComplexRegions.AbstractRegion):
                self.julia = julia_obj
        else:
            raise ValueError("Invalid argument to Region constructor")
        
    def get(self, field):
        return jl.getproperty(self.julia, jl.Symbol(field))

    def contains(self, z=None):
        if z is not None:
            return getattr(jl.ComplexRegions, "in")(z, self.julia)
        else:
            getattr(jl.ComplexRegions, "in")(self.julia)

    def boundary(self):
        b = jl.ComplexRegions.boundary(self.julia)
        return JuliaPath(b)
    
    def innerboundary(self):
        b = jl.ComplexRegions.innerboundary(self.julia)
        if isinstance(b, juliacall.VectorValue):  # type: ignore
            return [JuliaPath(j) for j in b]
        else:
            return JuliaPath(b)

    def outerboundary(self):
        b = jl.ComplexRegions.outerboundary(self.julia)
        if isinstance(b, juliacall.VectorValue):  # type: ignore
            paths = []
            for j in b:
                paths.append(JuliaPath(j))
            return paths
        else:
            return JuliaPath(b)

    def union(self, other):
        r = jl.ComplexRegions.union(self.julia, other.julia)
        return JuliaRegion(r)
    
    def intersect(self, other):
        r = jl.ComplexRegions.intersect(self.julia, other.julia)
        return JuliaRegion(r)
    
class ExteriorRegion(JuliaRegion):
    def __init__(self, inner):
        if isinstance(inner, juliacall.AnyValue):  # type: ignore
            if jl.isa(inner, jl.ComplexRegions.ExteriorRegion):
                self.julia = inner
            else:
                raise ValueError("Invalid argument to ExteriorRegion constructor")
        else:
            self.julia = jl.ComplexRegions.ExteriorRegion(inner)
        b = JuliaRegion.get(self, "inner")
        self.inner = [ClosedPath(j) for j in b]

    def isfinite(self):
        return False
    
    def __repr__(self):
        return f"Exterior region with {len(self.inner)} inner boundaries"
    
class InteriorRegion(JuliaRegion):
    def __init__(self, inner):
        if isinstance(inner, juliacall.AnyValue):  # type: ignore
            if jl.isa(inner, jl.ComplexRegions.InteriorConnectedRegion):
                self.julia = inner
            else:
                raise ValueError("Invalid argument to ExteriorRegion constructor")
        else:
            self.julia = jl.ComplexRegions.ExteriorRegion(inner)
        b = JuliaRegion.get(self, "inner")
        self.inner = [ClosedPath(j) for j in b]
        self.outer = JuliaRegion.get(self, "outer")

    def isfinite(self):
        return self.outer.isfinite() & all([b.isfinite() for b in self.inner])
    
    def __repr__(self):
        N = len(self.inner)
        if N==0:
            return f"Interior simply connected region"
        else:
            return f"Interior {N+1}-connnected region"
        
def between(curve1, curve2):
    """Construct the region between two closed curves."""
    r = jl.ComplexRegions.between(curve1.julia, curve2.julia)
    return InteriorRegion(r)

class Annulus(InteriorRegion):
    def __init__(self, outer, inner, center=0j):
        if isinstance(outer, juliacall.AnyValue):  # type: ignore
            if jl.isa(outer, jl.ComplexRegions.Annulus):
                self.julia = outer
            elif jl.isa(outer, jl.ComplexRegions.Circle) and jl.isa(inner, jl.ComplexRegions.Circle):
                self.julia = jl.ComplexRegions.Annulus(outer, inner)
            else:
                raise ValueError("Invalid argument to Annulus constructor")
        elif isinstance(inner, Circle) and isinstance(outer, Circle):
            self.julia = jl.ComplexRegions.Annulus(outer.julia, inner.julia)
        else:
            self.julia = jl.ComplexRegions.Annulus(outer, inner, center)

        self.inner = Circle(JuliaRegion.get(self, "inner"))
        self.outer = Circle(JuliaRegion.get(self, "outer"))

    def modulus(self):
        return jl.ComplexRegions.modulus(self.julia)
    
    def isfinite(self):
        return True

    def __repr__(self):
        return f"Annulus centered at {self.inner.center} with radii {self.inner.radius} and {self.outer.radius}"