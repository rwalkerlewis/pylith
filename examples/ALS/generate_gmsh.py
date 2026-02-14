#!/usr/bin/env nemesis
"""Generate a 2D mesh for the ALS example using Gmsh.

The mesh covers a rectangular domain based on the original ALS mesh coordinates.
We use the `gmsh_utils` module provided with PyLith for helper functions.

Run `generate_gmsh.py --help` to see the command line options.
Run `generate_gmsh.py --write` to generate the mesh.
"""

import gmsh
from pylith.meshio.gmsh_utils import (BoundaryGroup, MaterialGroup, GenerateMesh)


class App(GenerateMesh):
    """Application to generate the ALS example mesh using Gmsh.

    Domain covers approximately 19 km x 21 km based on original ALS data.
    Coordinates are in a projected coordinate system (e.g., UTM).

    p4-----------------p3
    |                   |
    |                   |
    |     material      |
    |                   |
    |                   |
    p1-----------------p2
    """
    # Domain boundaries based on original ALS mesh coordinates
    DOMAIN_W = 370000.0  # West boundary (m)
    DOMAIN_E = 389000.0  # East boundary (m)
    DOMAIN_S = 6527000.0  # South boundary (m)
    DOMAIN_N = 6548000.0  # North boundary (m)

    # Mesh discretization size (meters)
    DX = 1000.0

    def __init__(self):
        self.cell_choices = {
            "required": True,
            "choices": ["tri", "quad"],
        }
        self.filename = "mesh_tri.msh"

    def create_geometry(self):
        """Create geometry for the domain.
        """
        # Domain corner points
        p1 = gmsh.model.geo.add_point(self.DOMAIN_W, self.DOMAIN_S, 0.0)
        p2 = gmsh.model.geo.add_point(self.DOMAIN_E, self.DOMAIN_S, 0.0)
        p3 = gmsh.model.geo.add_point(self.DOMAIN_E, self.DOMAIN_N, 0.0)
        p4 = gmsh.model.geo.add_point(self.DOMAIN_W, self.DOMAIN_N, 0.0)

        # Domain boundary curves
        self.l_south = gmsh.model.geo.add_line(p1, p2)
        self.l_east = gmsh.model.geo.add_line(p2, p3)
        self.l_north = gmsh.model.geo.add_line(p3, p4)
        self.l_west = gmsh.model.geo.add_line(p4, p1)

        # Domain surface
        loop = gmsh.model.geo.add_curve_loop([self.l_south, self.l_east, self.l_north, self.l_west])
        self.s_domain = gmsh.model.geo.add_plane_surface([loop])

        gmsh.model.geo.synchronize()

    def mark(self):
        """Mark geometry for materials and boundary conditions.
        """
        # Define material
        materials = (
            MaterialGroup(tag=1, entities=[self.s_domain]),
        )
        for material in materials:
            material.create_physical_group()

        # Define boundary condition groups
        boundary_groups = (
            BoundaryGroup(name="boundary_xneg", tag=10, dim=1, entities=[self.l_west]),
            BoundaryGroup(name="boundary_xpos", tag=11, dim=1, entities=[self.l_east]),
            BoundaryGroup(name="boundary_yneg", tag=12, dim=1, entities=[self.l_south]),
            BoundaryGroup(name="boundary_ypos", tag=13, dim=1, entities=[self.l_north]),
        )
        for group in boundary_groups:
            group.create_physical_group()

    def generate_mesh(self, cell):
        """Generate the mesh.

        Args:
            cell: Cell type ("tri" or "quad")
        """
        gmsh.option.setNumber("Mesh.MeshSizeMin", self.DX)
        gmsh.option.setNumber("Mesh.MeshSizeMax", self.DX)

        if cell == "quad":
            gmsh.model.mesh.set_transfinite_automatic(recombine=True)
            self.filename = "mesh_quad.msh"
        else:
            gmsh.option.setNumber("Mesh.Algorithm", 8)
            self.filename = "mesh_tri.msh"

        gmsh.model.mesh.generate(2)
        gmsh.model.mesh.optimize("Laplace2D")


if __name__ == "__main__":
    App().main()


# End of file
