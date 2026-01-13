# cxregions - AI Development Context

## Project Overview

**cxregions** is a Python interface to the [ComplexRegions.jl](https://github.com/complexvariables/ComplexRegions.jl) Julia package, providing tools for working with complex regions defined by paths and curves. It enables users to create, manipulate, and analyze geometric shapes such as circular polygons and general polygons using Julia's computational capabilities from within Python.

### Key Information
- **Author**: Toby Driscoll (driscoll@udel.edu)
- **Version**: 0.1.2
- **License**: MIT
- **Python Requirements**: >= 3.11
- **Dependencies**: juliacall>=0.9.30,<0.10, numpy>=2.1,<3
- **Repository**: https://github.com/complexvariables/cxregions

## Architecture

The package uses [`juliacall`](https://github.com/JuliaPy/PythonCall.jl) to interface with Julia's ComplexRegions.jl package. The architecture follows a wrapper pattern where Python classes wrap Julia objects and provide Pythonic interfaces to Julia methods.

### Core Design Patterns

1. **Julia Object Wrapping**: Base classes like [`JuliaCurve`](src/cxregions/__init__.py:47) and [`JuliaPath`](src/cxregions/__init__.py:367) wrap Julia objects
2. Computation should be performed in Julia for efficiency as much as possible, while exposing a Pythonic API
3. **Type Hierarchy**: Specific curve/path types inherit from base wrappers
4. **Method Delegation**: Python methods delegate to corresponding Julia functions
5. **Automatic Type Conversion**: Complex numbers and arrays are automatically converted between Python and Julia
6. **Numpy Integration**: Where reasonable, results are returned as NumPy arrays for easy manipulation in Python

## Main Classes and Hierarchy

### Curves
- **Base**: [`JuliaCurve`](src/cxregions/__init__.py:47) → [`Curve`](src/cxregions/__init__.py:172)
- **Closed**: [`ClosedCurve`](src/cxregions/__init__.py:189) (adds winding number methods)
- **Specific Types**:
  - [`Line`](src/cxregions/__init__.py:205): Infinite straight lines
  - [`Segment`](src/cxregions/__init__.py:277): Line segments with endpoints
  - [`Circle`](src/cxregions/__init__.py:244): Circular curves
  - [`Arc`](src/cxregions/__init__.py:317): Circular arcs
  - [`Ray`](src/cxregions/__init__.py:302): Semi-infinite rays

### Paths
- **Base**: [`JuliaPath`](src/cxregions/__init__.py:367) → [`Path`](src/cxregions/__init__.py:513)
- **Closed**: [`ClosedPath`](src/cxregions/__init__.py:529) (adds containment methods)
- **Specific Types**:
  - [`Polygon`](src/cxregions/__init__.py:594): Polygons with straight sides
  - [`CircularPolygon`](src/cxregions/__init__.py:573): Polygons with circular arc sides
  - [`Rectangle`](src/cxregions/__init__.py:615): Axis-aligned rectangles

### Regions
- **Base**: [`JuliaRegion`](src/cxregions/__init__.py:664)
- **Types**:
  - [`Interior1CRegion`](src/cxregions/__init__.py:744): Simply connected interior regions
  - [`Exterior1CRegion`](src/cxregions/__init__.py:706): Simply connected exterior regions
  - [`InteriorConnectedRegion`](src/cxregions/__init__.py:762): Multiply connected interior regions
  - [`ExteriorRegion`](src/cxregions/__init__.py:724): Exterior regions with multiple boundaries
  - [`Annulus`](src/cxregions/__init__.py:819): Ring-shaped regions between two circles

## Common Methods and Properties

### Universal Methods (All Geometric Objects)
- [`point(t)`](src/cxregions/__init__.py:58): Get point at parameter t
- [`arclength()`](src/cxregions/__init__.py:62): Total arc length
- [`tangent(t)`](src/cxregions/__init__.py:65): Tangent vector at parameter t
- [`unittangent(t)`](src/cxregions/__init__.py:69): Unit tangent vector
- [`normal(t)`](src/cxregions/__init__.py:73): Normal vector
- [`conj()`](src/cxregions/__init__.py:80): Complex conjugate
- [`reverse()`](src/cxregions/__init__.py:84): Reverse orientation
- [`isfinite()`](src/cxregions/__init__.py:88): Check if finite
- [`dist(z)`](src/cxregions/__init__.py:117): Distance to point z
- [`closest(z)`](src/cxregions/__init__.py:114): Closest point to z
- [`reflect(z)`](src/cxregions/__init__.py:111): Reflect point z across curve

### Closed Curves/Paths Additional Methods
- [`winding(z)`](src/cxregions/__init__.py:199): Winding number around point z
- [`isinside(z)`](src/cxregions/__init__.py:546): Check if point is inside (ClosedPath only)

### Arithmetic Operations
All geometric objects support:
- Addition: `curve + complex_number`
- Subtraction: `curve - complex_number`
- Multiplication: `curve * complex_number`
- Division: `curve / complex_number`

## Testing Patterns and Examples

### Test Structure
Tests use pytest with fixtures for common objects. Key patterns:

#### 1. Fixture-Based Setup
```python
@pytest.fixture
def unit_circle():
    """Unit circle centered at origin."""
    return Circle(0, 1)

@pytest.fixture
def polygon():
    return Polygon([4, 4 + 3j, 3j, -2j, 6 - 2j, 6])
```

#### 2. Numerical Accuracy Testing
```python
def test_arclength(circle_thru_origin):
    c = circle_thru_origin
    assert abs(c.arclength() - 2 * np.sqrt(2) * np.pi) < 1e-7

def test_distance(circle_thru_origin):
    c = circle_thru_origin
    assert abs(c.dist(-1 + 1j) - np.sqrt(2)) < 1e-7
```

#### 3. Property Testing
```python
def test_line_not_finite(line_fixture):
    """Test that line is not finite"""
    assert not line_fixture.isfinite()

def test_line_is_positive(line_fixture):
    """Test that line is positive"""
    assert line_fixture.ispositive()
```

#### 4. Geometric Relationship Testing
```python
def test_polygon_winding(polygon):
    assert polygon.winding(5 - 1j) == 1
    assert polygon.winding(-1) == 0
    q = polygon.reverse()
    assert q.winding(5 - 1j) == -1
```

#### 5. Arithmetic Operation Testing
```python
def test_circle_add(circle_thru_origin):
    l1 = circle_thru_origin
    l2 = l1 + 2j
    assert l2.dist(l1.point(0.5) + 2j) < 1e-7
    l2 = 4 - 1j + l1
    assert l2.dist(l1.point(0.5) + (4 - 1j)) < 1e-7
```

#### 6. Type and Instance Testing
```python
def test_line_inverse_is_circle(line_fixture):
    """Test that 1/line is a Circle"""
    inv_l = line_fixture.inv()
    assert isinstance(inv_l, Circle)

def test_annulus_from_circles(small_circle, large_circle):
    """Test creating annulus from two circles."""
    annulus = Annulus(large_circle, small_circle)
    assert isinstance(annulus, Annulus)
    assert annulus.isfinite()
```

#### 7. Intersection Testing
```python
def test_intersect_circles():
    u = 1/5 + 1j/2
    c = Circle(0, 1)
    z = c.intersect(Circle(u, 3/2))
    assert isinstance(z, np.ndarray)
    assert len(z) == 2
    assert np.allclose(np.abs(z - 0), 1)
```

#### 8. Region Containment Testing
```python
def test_annulus_containment(small_circle, large_circle):
    """Test point containment in annulus."""
    annulus = Annulus(large_circle, small_circle)
    
    # Point in annulus (between circles)
    assert annulus.contains(1)
    
    # Point inside inner circle
    assert not annulus.contains(0.25)
    
    # Point outside outer circle
    assert not annulus.contains(3)
```

### Helper Functions in Tests
```python
def check(u, v, rtol=1e-12, atol=1e-12):
    """Helper function to check approximate equality with tolerance"""
    if isinstance(u, (list, tuple, np.ndarray)) and isinstance(v, (list, tuple, np.ndarray)):
        return all(check(ui, vi, rtol, atol) for ui, vi in zip(u, v))
    elif isinstance(u, (list, tuple, np.ndarray)):
        return all(check(ui, v, rtol, atol) for ui in u)
    else:
        return np.isclose(u, v, rtol=rtol, atol=atol)
```

## Common Usage Patterns

### Creating Basic Shapes
```python
from cxregions import Circle, Line, Polygon, Arc, Segment

# Circle: center and radius
circle = Circle(0, 1)  # unit circle at origin
circle = Circle(1+1j, 2)  # circle at (1,1) with radius 2

# Line: two points or point + direction
line = Line(0, 1+1j)  # line through origin and (1,1)
line = Line(0, direction=1+1j)  # line through origin in direction (1,1)

# Polygon from vertices
poly = Polygon([0, 1, 1+1j, 1j])  # unit square

# Arc: start, end, center
arc = Arc(1, 1j, 0)  # quarter circle from 1 to i around origin
```

### Working with Regions
```python
from cxregions import interior, exterior, disk, Annulus

# Create regions
inside_circle = interior(Circle(0, 1))
outside_circle = exterior(Circle(0, 1))
disk_region = disk(0, 1)  # same as interior(Circle(0, 1))
ring = Annulus(2, 1, center=0)  # annulus with outer radius 2, inner radius 1

# Test containment
assert inside_circle.contains(0.5)
assert not inside_circle.contains(2)
assert ring.contains(1.5)  # between inner and outer circles
```

### Geometric Operations
```python
# Point evaluation
z = circle.point(0.25)  # point at parameter t=0.25

# Geometric properties
length = curve.arclength()
tangent_vec = curve.tangent(0.5)
distance = curve.dist(1+1j)
closest_pt = curve.closest(1+1j)

# Transformations
translated = curve + 1+1j  # translate by (1,1)
scaled = curve * 2  # scale by factor 2
rotated = curve * 1j  # rotate by 90 degrees
```

## Development Guidelines

### When Writing Tests
1. **Use fixtures** for commonly used geometric objects
2. **Test numerical accuracy** with appropriate tolerances (typically 1e-7 to 1e-12)
3. **Test both positive and negative cases** for containment and properties
4. **Verify type consistency** after operations
5. **Test edge cases** like boundary points, degenerate cases
6. **Use descriptive test names** that explain what is being tested
7. **Group related tests** in classes (e.g., `TestAnnulus`, `TestRegionBoundaries`)

### Common Test Assertions
- `assert isinstance(obj, ExpectedType)` for type checking
- `assert np.isclose(actual, expected, rtol=1e-7)` for floating point comparisons
- `assert np.allclose(array1, array2)` for array comparisons
- `assert obj.contains(point)` or `assert not obj.contains(point)` for containment
- `assert obj.isfinite()` or `assert not obj.isfinite()` for finiteness

### Error Handling
The package generally raises `ValueError` for invalid constructor arguments or incompatible operations. Tests should verify both successful operations and appropriate error conditions.

## File Structure
```
cxregions/
├── src/cxregions/__init__.py    # Main package implementation
├── tests/                       # Test suite
│   ├── test_regions.py         # Region and containment tests
│   ├── test_circle.py          # Circle-specific tests
│   ├── test_line.py            # Line-specific tests
│   ├── test_polygon.py         # Polygon tests
│   ├── test_intersections.py   # Intersection tests
│   └── ...                     # Other curve/path tests
├── pyproject.toml              # Project configuration
├── README.md                   # User documentation
└── setup.py                    # Legacy setup file
```

This context should provide sufficient information for an AI to understand the codebase structure, testing patterns, and development practices for the cxregions package.