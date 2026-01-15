"""
Möbius transformation class for the cxregions package.
"""

import juliacall
import numpy as np
from . import jl, JLCR
from .curves import wrap_jl_curve
from .paths import wrap_jl_path
from .regions import wrap_jl_region

class Mobius:
    """
    Representation of a Möbius or bilinear transformation.
    
    A Möbius transformation is a function of the form f(z) = (az + b) / (cz + d).
    """
    
    def __init__(self, *args):
        """
        Initialize a Möbius transformation.
        
        Parameters
        ----------
        *args : 
            - Mobius(a, b, c, d): Construct from coefficients
            - Mobius(A): Construct from 2x2 matrix
            - Mobius(source, image): Construct map from 3 points to 3 points
            - Mobius(c1, c2): Construct map from one curve/line to another
            - Mobius(julia_obj): Wrap an existing Julia Mobius object

        Raises
        ------
        ValueError
            If the arguments do not match any valid constructor.

        Examples
        --------
        >>> f = Mobius(2, 1, 1, 2)  # f(z) = (2z + 1) / (z + 2)
        >>> f = Mobius([[2, 1], [1, 2]])
        >>> f = Mobius([0, 1, np.inf], [1, 1j, -1])  # Map 0,1,inf to 1,i,-1
        >>> f = Mobius(Line(-1, 1), Circle(0, 1))
        """
        if len(args) == 1 and isinstance(args[0], juliacall.AnyValue): # type: ignore
            if jl.isa(args[0], JLCR.Mobius):
                self.julia = args[0]
                return
        
        # Helper to get julia objects from potential wrappers
        def get_jl(x):
            return getattr(x, 'julia', x)

        try:
            if len(args) == 4:
                self.julia = JLCR.Mobius(*args)
            elif len(args) == 2:
                self.julia = JLCR.Mobius(get_jl(args[0]), get_jl(args[1]))
            elif len(args) == 1:
                # Could be a matrix (numpy array or list of lists)
                self.julia = JLCR.Mobius(np.asarray(args[0]))
            else:
                raise ValueError("Invalid number of arguments for Mobius constructor")
        except Exception as e:
            raise ValueError(f"Failed to construct Mobius: {e}")

    @property
    def coeff(self):
        """
        Return the coefficients [a, b, c, d] of the Möbius transformation.
        """
        return np.array(self.julia.coeff)

    def __call__(self, z):
        """
        Evaluate the Möbius transformation.
        
        Parameters
        ----------
        z : complex, array_like, Curve, Path, or Region
            Input to the transformation
            
        Returns
        -------
        Same type as z
            The image of z under the transformation
        """
        if isinstance(z, (list, tuple, np.ndarray)):
            # Use broadcasting in Julia for arrays of points
            res = jl.broadcast(self.julia, z)
            return np.array(res)
        
        j_z = getattr(z, 'julia', z)
        res = self.julia(j_z)
        
        if jl.isa(res, JLCR.AbstractCurve):
            return wrap_jl_curve(res)
        elif jl.isa(res, JLCR.AbstractPath):
            return wrap_jl_path(res)
        elif jl.isa(res, JLCR.AbstractRegion):
            return wrap_jl_region(res)
        else:
            # For complex numbers, juliacall should return Python types
            return res

    def inv(self):
        """
        Return the inverse Möbius transformation.
        """
        return Mobius(jl.inv(self.julia))
    
    def __matmul__(self, other):
        """
        Compose two Möbius transformations using the @ operator.
        """
        if isinstance(other, Mobius):
            # In Julia, composition is ∘
            compose = getattr(jl, "∘")
            return Mobius(compose(self.julia, other.julia))
        return NotImplemented

    def compose(self, other):
        """
        Compose two Möbius transformations.
        """
        if isinstance(other, Mobius):
            return Mobius(getattr(jl, "∘")(self.julia, other.julia))
        return NotImplemented

    def __repr__(self):
        return f"Mobius({jl.repr(self.julia)})"

    def __str__(self):
        # Clean up the output slightly for strings if needed, 
        # but Julia's show for Mobius is nice.
        return str(self.julia)
