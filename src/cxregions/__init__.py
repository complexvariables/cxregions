"""
cxregions: A Python interface to ComplexRegions.jl

This package provides tools for working with complex regions defined by paths and curves.
It allows users to create, manipulate, and analyze geometric shapes such as circular 
polygons and general polygons using Julia's computational capabilities from within Python.
"""

import juliacall
import numpy as np

# Initialize Julia environment and load ComplexRegions.jl
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

# Julia ComplexRegions module reference
JLCR = jl.ComplexRegions

# Import all classes and functions from submodules
from .curves import (
    JuliaCurve, Curve, ClosedCurve, Line, Circle, Segment, Ray, Arc,
    wrap_jl_curve, unitcircle
)

from .mobius import Mobius

from .paths import (
    JuliaPath, Path, ClosedPath, CircularPolygon, Polygon, Rectangle,
    wrap_jl_path, n_gon, quad
)

from .regions import (
    JuliaRegion, Exterior1CRegion, Interior1CRegion, ExteriorRegion, 
    InteriorConnectedRegion, Annulus, wrap_jl_region, Jordan, get_julia,
    between, interior, exterior, disk, unitdisk, halfplane,
    upperhalfplane, lowerhalfplane, lefthalfplane, righthalfplane
)

# Define public API
__all__ = [
    # Julia interface
    "jl", "update",
    
    # Curve classes
    "Curve", "ClosedCurve", "Line", "Segment", "Circle", "Ray", "Arc",
    
    # Path classes  
    "Path", "ClosedPath", "CircularPolygon", "Polygon", "Rectangle",
    
    # Transformation classes
    "Mobius",

    # Region classes
    "Exterior1CRegion", "Interior1CRegion", "ExteriorRegion", "InteriorConnectedRegion",
    "Annulus",
    
    # Utilities
    "n_gon", "unitcircle", "unitdisk", "between", "interior", "exterior", "disk", "quad",
    "halfplane", "upperhalfplane", "lowerhalfplane", "lefthalfplane", "righthalfplane"
]

def update():
    """Check for and install updates to ComplexRegions.jl.

    Changes take effect on the next Python session.
    """
    before = {v.name: str(v.version) for v in jl.Pkg.dependencies().values()
              if v.name == "ComplexRegions"}
    jl.seval('Pkg.update("ComplexRegions")')
    after = {v.name: str(v.version) for v in jl.Pkg.dependencies().values()
             if v.name == "ComplexRegions"}
    old = before.get("ComplexRegions", "unknown")
    new = after.get("ComplexRegions", "unknown")
    if old == new:
        print(f"ComplexRegions.jl is already up to date (v{new}).")
    else:
        print(f"ComplexRegions.jl updated from v{old} to v{new}. Restart Python to apply.")


# Package metadata
from importlib.metadata import version as _version, PackageNotFoundError as _PackageNotFoundError
try:
    __version__ = _version("cxregions")
except _PackageNotFoundError:
    __version__ = "unknown"
__author__ = "Toby Driscoll"
__email__ = "driscoll@udel.edu"
__description__ = "A Python interface to the ComplexRegions.jl Julia package"