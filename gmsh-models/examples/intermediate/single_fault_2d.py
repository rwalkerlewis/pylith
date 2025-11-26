#!/usr/bin/env nemesis
"""Generate a 2D mesh with a single vertical strike-slip fault using Gmsh.

This intermediate example demonstrates:
- Domain with a fault that splits the domain
- Two materials (one on each side of the fault)
- Fault marking for PyLith
- Variable mesh refinement near the fault

Domain: 100km x 150km
-50.0 km <= x <= 50.0 km
-75.0 km <= y <= 75.0 km

Fault: Vertical line along x=0 (y-axis)

Run `single_fault_2d.py --write` to generate the mesh.
"""

import gmsh
from pylith.meshio.gmsh_utils import (BoundaryGroup, MaterialGroup, GenerateMesh)


class App(GenerateMesh):
    """
    Application used to generate a 2D mesh with a single vertical fault.
    
    This example shows:
    - Domain split by a fault
    - Multiple materials
    - Fault marking
    - Mesh refinement near fault
    """
    
    DOMAIN_X = 100.0e+3  # 100 km
    DOMAIN_Y = 150.0e+3  # 150 km
    
    DX_FAULT = 5.0e+3    # Cell size on fault: 5 km
    DX_BIAS = 1.15       # Bias factor for cell size growth
    
    def __init__(self):
        """Constructor.
        """
        super().__init__()
        self.cell_choices = {
            "default": "tri",
            "choices": ["tri", "quad"],
        }
        self.filename = "single_fault_2d.msh"
    
    def create_geometry(self):
        """Create geometry with a vertical fault splitting the domain.
        """
        lx = self.DOMAIN_X
        ly = self.DOMAIN_Y
        x1 = -0.5 * lx
        y1 = -0.5 * ly
        
        # Create corner points
        p1 = gmsh.model.geo.add_point(x1, y1, 0.0)
        p2 = gmsh.model.geo.add_point(x1 + lx, y1, 0.0)
        p3 = gmsh.model.geo.add_point(x1 + lx, y1 + ly, 0.0)
        p4 = gmsh.model.geo.add_point(x1, y1 + ly, 0.0)
        
        # Create points on fault (at x=0)
        p5 = gmsh.model.geo.add_point(0.0, y1, 0.0)
        p6 = gmsh.model.geo.add_point(0.0, y1 + ly, 0.0)
        
        # Create curves
        self.c_yneg1 = gmsh.model.geo.add_line(p1, p5)
        self.c_yneg2 = gmsh.model.geo.add_line(p5, p2)
        self.c_xpos = gmsh.model.geo.add_line(p2, p3)
        self.c_ypos2 = gmsh.model.geo.add_line(p3, p6)
        self.c_ypos1 = gmsh.model.geo.add_line(p6, p4)
        self.c_xneg = gmsh.model.geo.add_line(p4, p1)
        self.c_fault = gmsh.model.geo.add_line(p5, p6)
        
        # Create surfaces (left and right of fault)
        loop0 = gmsh.model.geo.add_curve_loop([self.c_yneg1, self.c_fault, self.c_ypos1, self.c_xneg])
        self.s_left = gmsh.model.geo.add_plane_surface([loop0])
        
        loop1 = gmsh.model.geo.add_curve_loop([self.c_yneg2, self.c_xpos, self.c_ypos2, -self.c_fault])
        self.s_right = gmsh.model.geo.add_plane_surface([loop1])
        
        gmsh.model.geo.synchronize()
    
    def mark(self):
        """Mark geometry for materials, boundary conditions, and fault.
        """
        # Create two materials (one on each side of fault)
        materials = (
            MaterialGroup(tag=1, entities=[self.s_left]),
            MaterialGroup(tag=2, entities=[self.s_right]),
        )
        for material in materials:
            material.create_physical_group()
        
        # Create boundary groups
        boundary_groups = (
            BoundaryGroup(name="boundary_xneg", tag=10, dim=1, entities=[self.c_xneg]),
            BoundaryGroup(name="boundary_xpos", tag=11, dim=1, entities=[self.c_xpos]),
            BoundaryGroup(name="boundary_yneg", tag=12, dim=1, entities=[self.c_yneg1, self.c_yneg2]),
            BoundaryGroup(name="boundary_ypos", tag=13, dim=1, entities=[self.c_ypos1, self.c_ypos2]),
            BoundaryGroup(name="fault", tag=20, dim=1, entities=[self.c_fault]),
        )
        for group in boundary_groups:
            group.create_physical_group()
    
    def generate_mesh(self, cell):
        """Generate mesh with refinement near the fault.
        """
        # Turn off default sizing methods
        gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)
        
        # Create distance field from fault
        field_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(field_distance, "CurvesList", [self.c_fault])
        
        # Create size field with geometric progression
        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(field_distance, min_dx=self.DX_FAULT, bias=self.DX_BIAS)
        gmsh.model.mesh.field.setString(field_size, "F", math_exp)
        
        # Use the size field as background mesh
        gmsh.model.mesh.field.setAsBackgroundMesh(field_size)
        
        if cell == "quad":
            gmsh.option.setNumber("Mesh.Algorithm", 8)
            gmsh.model.mesh.generate(2)
            gmsh.model.mesh.recombine()
        else:
            gmsh.model.mesh.generate(2)
        
        gmsh.model.mesh.optimize("Laplace2D")


if __name__ == "__main__":
    App().main()
