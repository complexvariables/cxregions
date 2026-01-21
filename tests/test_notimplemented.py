"""
Tests for proper NotImplemented behavior in magic methods.

This module tests that arithmetic magic methods return NotImplemented
for unsupported types, allowing Python to raise appropriate TypeErrors
instead of propagating Julia exceptions.
"""

import pytest
from cxregions import Circle, Line, Segment, Arc, Ray, Polygon, Rectangle, CircularPolygon


class TestCurveNotImplemented:
    """Test that curve arithmetic methods return NotImplemented for unsupported types."""
    
    def test_circle_add_unsupported(self):
        """Test Circle + unsupported type raises TypeError."""
        c = Circle(0, 1)
        with pytest.raises(TypeError, match="unsupported operand type"):
            c + "invalid"
    
    def test_circle_radd_unsupported(self):
        """Test unsupported type + Circle raises TypeError."""
        c = Circle(0, 1)
        with pytest.raises(TypeError):
            "invalid" + c
    
    def test_circle_sub_unsupported(self):
        """Test Circle - unsupported type raises TypeError."""
        c = Circle(0, 1)
        with pytest.raises(TypeError, match="unsupported operand type"):
            c - "invalid"
    
    def test_circle_rsub_unsupported(self):
        """Test unsupported type - Circle raises TypeError."""
        c = Circle(0, 1)
        with pytest.raises(TypeError):
            "invalid" - c
    
    def test_circle_mul_unsupported(self):
        """Test Circle * unsupported type raises TypeError."""
        c = Circle(0, 1)
        with pytest.raises(TypeError):
            c * "invalid"
    
    def test_circle_rmul_unsupported(self):
        """Test unsupported type * Circle raises TypeError."""
        c = Circle(0, 1)
        with pytest.raises(TypeError):
            "invalid" * c
    
    def test_circle_truediv_unsupported(self):
        """Test Circle / unsupported type raises TypeError."""
        c = Circle(0, 1)
        with pytest.raises(TypeError):
            c / "invalid"
    
    def test_line_add_unsupported(self):
        """Test Line + unsupported type raises TypeError."""
        line = Line(0, 1)
        with pytest.raises(TypeError, match="unsupported operand type"):
            line + "invalid"
    
    def test_segment_mul_unsupported(self):
        """Test Segment * unsupported type raises TypeError."""
        seg = Segment(0, 1)
        with pytest.raises(TypeError):
            seg * "invalid"
    
    def test_arc_div_unsupported(self):
        """Test Arc / unsupported type raises TypeError."""
        arc = Arc(1, 1j, 0)
        with pytest.raises(TypeError):
            arc / "invalid"
    
    def test_ray_sub_unsupported(self):
        """Test Ray - unsupported type raises TypeError."""
        import numpy as np
        ray = Ray(0, np.pi/4)
        with pytest.raises(TypeError, match="unsupported operand type"):
            ray - "invalid"


class TestCurveValidOperations:
    """Test that valid arithmetic operations still work correctly."""
    
    def test_circle_add_complex(self):
        """Test Circle + complex number works."""
        c = Circle(0, 1)
        c2 = c + (1 + 1j)
        assert c2.center == (1 + 1j)
        assert c2.radius == 1
    
    def test_circle_radd_complex(self):
        """Test complex number + Circle works."""
        c = Circle(0, 1)
        c2 = (1 + 1j) + c
        assert c2.center == (1 + 1j)
        assert c2.radius == 1
    
    def test_circle_mul_scalar(self):
        """Test Circle * scalar works."""
        c = Circle(0, 1)
        c2 = c * 2
        assert c2.radius == 2
    
    def test_circle_rmul_scalar(self):
        """Test scalar * Circle works."""
        c = Circle(0, 1)
        c2 = 2 * c
        assert c2.radius == 2
    
    def test_circle_div_scalar(self):
        """Test Circle / scalar works."""
        c = Circle(0, 2)
        c2 = c / 2
        assert c2.radius == 1
    
    def test_line_add_complex(self):
        """Test Line + complex number works."""
        line = Line(0, direction=1)
        line2 = line + (1 + 1j)
        assert abs(line2.base - (1 + 1j)) < 1e-10
    
    def test_segment_mul_complex(self):
        """Test Segment * complex number works."""
        seg = Segment(0, 1)
        seg2 = seg * 2
        assert abs(seg2.last - 2) < 1e-10


class TestPathNotImplemented:
    """Test that path arithmetic methods return NotImplemented for unsupported types."""
    
    def test_polygon_add_unsupported(self):
        """Test Polygon + unsupported type raises TypeError."""
        p = Polygon([0, 1, 1 + 1j, 1j])
        with pytest.raises(TypeError, match="unsupported operand type"):
            p + "invalid"
    
    def test_polygon_radd_unsupported(self):
        """Test unsupported type + Polygon raises TypeError."""
        p = Polygon([0, 1, 1 + 1j, 1j])
        with pytest.raises(TypeError):
            "invalid" + p
    
    def test_polygon_sub_unsupported(self):
        """Test Polygon - unsupported type raises TypeError."""
        p = Polygon([0, 1, 1 + 1j, 1j])
        with pytest.raises(TypeError, match="unsupported operand type"):
            p - "invalid"
    
    def test_polygon_rsub_unsupported(self):
        """Test unsupported type - Polygon raises TypeError."""
        p = Polygon([0, 1, 1 + 1j, 1j])
        with pytest.raises(TypeError):
            "invalid" - p
    
    def test_polygon_mul_unsupported(self):
        """Test Polygon * unsupported type raises TypeError."""
        p = Polygon([0, 1, 1 + 1j, 1j])
        with pytest.raises(TypeError):
            p * "invalid"
    
    def test_polygon_rmul_unsupported(self):
        """Test unsupported type * Polygon raises TypeError."""
        p = Polygon([0, 1, 1 + 1j, 1j])
        with pytest.raises(TypeError):
            "invalid" * p
    
    def test_polygon_truediv_unsupported(self):
        """Test Polygon / unsupported type raises TypeError."""
        p = Polygon([0, 1, 1 + 1j, 1j])
        with pytest.raises(TypeError):
            p / "invalid"
    
    def test_rectangle_add_unsupported(self):
        """Test Rectangle + unsupported type raises TypeError."""
        import numpy as np
        rect = Rectangle(0 + 0j, np.array([1.0, 0.5]))
        with pytest.raises(TypeError, match="unsupported operand type"):
            rect + "invalid"
    
    def test_circularpolygon_mul_unsupported(self):
        """Test CircularPolygon * unsupported type raises TypeError."""
        from cxregions import Segment
        seg1 = Segment(0, 1)
        seg2 = Segment(1, 1 + 1j)
        seg3 = Segment(1 + 1j, 1j)
        seg4 = Segment(1j, 0)
        cpoly = CircularPolygon([seg1, seg2, seg3, seg4])
        with pytest.raises(TypeError):
            cpoly * "invalid"


class TestPathValidOperations:
    """Test that valid path arithmetic operations still work correctly."""
    
    def test_polygon_add_complex(self):
        """Test Polygon + complex number works."""
        p = Polygon([0, 1, 1 + 1j, 1j])
        p2 = p + (1 + 1j)
        assert abs(p2.vertex(1) - (1 + 1j)) < 1e-10
    
    def test_polygon_radd_complex(self):
        """Test complex number + Polygon works."""
        p = Polygon([0, 1, 1 + 1j, 1j])
        p2 = (1 + 1j) + p
        assert abs(p2.vertex(1) - (1 + 1j)) < 1e-10
    
    def test_polygon_mul_scalar(self):
        """Test Polygon * scalar works."""
        p = Polygon([0, 1, 1 + 1j, 1j])
        p2 = p * 2
        assert abs(p2.vertex(2) - 2) < 1e-10
    
    def test_polygon_rmul_scalar(self):
        """Test scalar * Polygon works."""
        p = Polygon([0, 1, 1 + 1j, 1j])
        p2 = 2 * p
        assert abs(p2.vertex(2) - 2) < 1e-10
    
    def test_polygon_div_scalar(self):
        """Test Polygon / scalar works."""
        p = Polygon([0, 2, 2 + 2j, 2j])
        p2 = p / 2
        assert abs(p2.vertex(2) - 1) < 1e-10
    
    def test_rectangle_add_complex(self):
        """Test Rectangle + complex number works."""
        import numpy as np
        rect = Rectangle(0 + 0j, np.array([1.0, 0.5]))
        rect2 = rect + (1 + 1j)
        assert abs(rect2.center - (1 + 1j)) < 1e-10


class TestEdgeCases:
    """Test edge cases and special scenarios."""
    
    def test_none_type(self):
        """Test that None raises TypeError."""
        c = Circle(0, 1)
        with pytest.raises(TypeError):
            c + None
    
    def test_list_type(self):
        """Test that list raises TypeError."""
        c = Circle(0, 1)
        with pytest.raises(TypeError):
            c + [1, 2, 3]
    
    def test_dict_type(self):
        """Test that dict raises TypeError."""
        c = Circle(0, 1)
        with pytest.raises(TypeError):
            c * {"key": "value"}
    
    def test_mixed_curve_types(self):
        """Test that operations between different curve types work if Julia supports them."""
        # This should work - adding a complex number to any curve
        c = Circle(0, 1)
        line = Line(0, 1)
        seg = Segment(0, 1)
        
        # All should accept complex addition
        c2 = c + 1j
        line2 = line + 1j
        seg2 = seg + 1j
        
        assert isinstance(c2, Circle)
        assert isinstance(line2, Line)
        assert isinstance(seg2, Segment)
