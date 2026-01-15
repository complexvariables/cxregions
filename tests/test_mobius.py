import pytest
import numpy as np
from cxregions import Mobius, Circle, Line, interior, disk, upperhalfplane

def test_mobius_coeffs():
    # f(z) = 2z + 1 / (z + 2)
    f = Mobius(2, 1, 1, 2)
    assert f(0) == pytest.approx(0.5)
    assert f(1) == pytest.approx(1.0)
    assert f(np.inf) == pytest.approx(2.0)

def test_mobius_3pt():
    # Map 0, 1, inf to 1, i, -1 (should be Cayley-ish or similar)
    f = Mobius([0, 1, np.inf], [1, 1j, -1])
    assert f(0) == pytest.approx(1)
    assert f(1) == pytest.approx(1j)
    assert f(np.inf) == pytest.approx(-1)
    assert f([0, 1, np.inf]) == pytest.approx([1, 1j, -1])

def test_mobius_curve():
    # Map upper half plane boundary (Line) to unit circle
    f = Mobius(Line(-1, 1), Circle(0, 1))
    c = f(Line(-1, 1))
    assert isinstance(c, Circle)
    assert c.center == pytest.approx(0)
    assert c.radius == pytest.approx(1)

def test_mobius_region():
    # Map upper half plane to unit disk
    # This might require a specific Mobius map
    f = Mobius([1j, 1, -1], [0, 1, 1j]) # 1j->0, 1->1, -1->1j
    d = f(upperhalfplane)
    assert d.contains(f(2j)) # type: ignore
    # upperhalfplane boundary is the real line. 
    # f(real line) should pass through f(0)=inf, so it's a Line
    b = f(Line(-1, 1))
    from cxregions import Line as RegLine
    assert isinstance(b, RegLine)

def test_mobius_inverse():
    f = Mobius(2, 3, 4, 5)
    finv = f.inv()
    z = 1 + 1j
    assert finv(f(z)) == pytest.approx(z)

def test_mobius_composition():
    f = Mobius(1, 2, 3, 4)
    g = Mobius(5, 6, 7, 8)
    h = f @ g
    z = 2 - 3j
    assert isinstance(h, Mobius)
    assert h(z) == pytest.approx(f(g(z))) # type: ignore

def test_mobius_array():
    f = Mobius(1, 0, 0, 1) # identity
    z = np.array([1, 2, 3])
    w = f(z)
    assert z == pytest.approx(w)
