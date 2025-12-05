#!/usr/bin/env nemesis
"""Generate a simple 2D box mesh with no faults using Gmsh.

This is the most basic example - a rectangular domain with boundary conditions
but no faults. This demonstrates the fundamental structure needed for PyLith
mesh generation.

Domain: 100km x 100km
-50.0 km <= x <= 50.0 km
-50.0 km <= y <= 50.0 km

Run `box_2d_no_faults.py --write` to generate the mesh.
"""

import gmsh
from pylith.meshio.gmsh_utils import (BoundaryGroup, MaterialGroup, GenerateMesh)


class App(GenerateMesh):
    """
    Application used to generate a simple 2D box mesh with no faults.
    
    This is the most basic example showing:
    - Domain creation
    - Material marking
    - Boundary condition marking
    - Mesh generation
    """
    
    DOMAIN_X = 100.0e+3  # 100 km
    DOMAIN_Y = 100.0e+3  # 100 km
    
    DX_DEFAULT = 10.0e+3  # Default cell size: 10 km
    
    def __init__(self):
        """Constructor.
        """
        super().__init__()
        self.cell_choices = {
            "default": "tri",
            "choices": ["tri", "quad"],
        }
        self.filename = "box_2d_no_faults.msh"
    
    def create_geometry(self):
        """Create geometry - a simple rectangular domain.
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
        
        # Create boundary curves
        self.c_xneg = gmsh.model.geo.add_line(p1, p4)
        self.c_xpos = gmsh.model.geo.add_line(p2, p3)
        self.c_yneg = gmsh.model.geo.add_line(p1, p2)
        self.c_ypos = gmsh.model.geo.add_line(p4, p3)
        
        # Create surface
        loop = gmsh.model.geo.add_curve_loop([self.c_xneg, self.c_ypos, -self.c_xpos, -self.c_yneg])
        self.s_domain = gmsh.model.geo.add_plane_surface([loop])
        
        gmsh.model.geo.synchronize()
    
    def mark(self):
        """Mark geometry for materials and boundary conditions.
        """
        # Create a single material for the entire domain
        materials = (
            MaterialGroup(tag=1, entities=[self.s_domain]),
        )
        for material in materials:
            material.create_physical_group()
        
        # Create boundary groups
        boundary_groups = (
            BoundaryGroup(name="boundary_xneg", tag=10, dim=1, entities=[self.c_xneg]),
            BoundaryGroup(name="boundary_xpos", tag=11, dim=1, entities=[self.c_xpos]),
            BoundaryGroup(name="boundary_yneg", tag=12, dim=1, entities=[self.c_yneg]),
            BoundaryGroup(name="boundary_ypos", tag=13, dim=1, entities=[self.c_ypos]),
        )
        for group in boundary_groups:
            group.create_physical_group()
    
    def generate_mesh(self, cell):
        """Generate the mesh with uniform cell size.
        """
        # Set uniform cell size
        gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.set_number("Mesh.MeshSizeMax", self.DX_DEFAULT)
        gmsh.option.set_number("Mesh.MeshSizeMin", self.DX_DEFAULT)
        
        if cell == "quad":
            gmsh.option.setNumber("Mesh.Algorithm", 8)
            gmsh.model.mesh.generate(2)
            gmsh.model.mesh.recombine()
        else:
            gmsh.model.mesh.generate(2)
        
        gmsh.model.mesh.optimize("Laplace2D")


if __name__ == "__main__":
    App().main()
