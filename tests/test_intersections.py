import pytest
import numpy as np
from cxregions import Line, Segment, Arc, Circle, Ray

def test_intersect_circles():
    u = 1/5 + 1j/2
    c = Circle(0, 1)
    z = c.intersect(Circle(u, 3/2))
    assert isinstance(z, np.ndarray)
    assert len(z) == 2
    assert np.allclose(np.abs(z - 0), 1)
    assert np.allclose(np.abs(z - u), 3/2)
    z = c.intersect(Circle(u, 0.1))
    assert isinstance(z, np.ndarray) and len(z) == 0

# TODO more tests for intersections