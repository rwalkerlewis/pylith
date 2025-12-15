#!/usr/bin/env nemesis
"""Generate a tri mesh for thermoelasticity examples using Gmsh.

This mesh represents a vertical cross-section through the crust with an
optional strike-slip fault. The domain is designed to demonstrate:
- Thermoelastic coupling
- Thermal stresses from temperature changes
- Heat conduction

Domain: 40 km wide x 20 km deep
- -20.0 km <= x <= 20.0 km
- -20.0 km <= y <= 0.0 km (y=0 is the surface)

The fault runs vertically through the center of the domain from the surface
down to 15 km depth (buried at the bottom).

        surface (y=0)
    p4-----p6-----p3
    |       |      |
    |       |(fault|
    |       |      |
    |      p7      |  <- fault tip at 15 km depth
    |              |
    p1------------p2
        bottom (y=-20 km)

Run `generate_gmsh.py --help` to see command line options.
Run `generate_gmsh.py --write` to generate the mesh.
"""

import gmsh

from pylith.meshio.gmsh_utils import (BoundaryGroup, MaterialGroup, GenerateMesh, VertexGroup)


class App(GenerateMesh):
    """Application for generating mesh with Gmsh.
    """
    DOMAIN_X = 40.0e+3  # 40 km wide
    DOMAIN_Y = 20.0e+3  # 20 km deep
    FAULT_DEPTH = 15.0e+3  # Fault extends to 15 km depth

    DX_FAULT = 500.0  # 500 m resolution near fault
    DX_BIAS = 1.15  # Geometric progression bias

    def __init__(self):
        """Constructor.
        """
        self.cell_choices = {
            "default": "tri",
            "choices": ["tri", "quad"],
        }
        self.filename = "mesh_tri.msh"

    def create_geometry(self):
        """Create geometry.
        """
        lx = self.DOMAIN_X
        ly = self.DOMAIN_Y
        fault_depth = self.FAULT_DEPTH

        x1 = -0.5 * lx  # Left edge
        y1 = -ly  # Bottom edge (depth)

        # Create corner points
        p1 = gmsh.model.geo.add_point(x1, y1, 0.0)  # Bottom-left
        p2 = gmsh.model.geo.add_point(x1 + lx, y1, 0.0)  # Bottom-right
        p3 = gmsh.model.geo.add_point(x1 + lx, 0.0, 0.0)  # Top-right (surface)
        p4 = gmsh.model.geo.add_point(x1, 0.0, 0.0)  # Top-left (surface)

        # Fault points (vertical fault at x=0)
        p5 = gmsh.model.geo.add_point(0.0, y1, 0.0)  # Fault at bottom (for mesh transition)
        p6 = gmsh.model.geo.add_point(0.0, 0.0, 0.0)  # Fault at surface
        p7 = gmsh.model.geo.add_point(0.0, -fault_depth, 0.0)  # Fault tip (buried edge)

        # Create curves
        # Bottom boundary (in two parts for mesh transition)
        self.c_bottom1 = gmsh.model.geo.add_line(p1, p5)
        self.c_bottom2 = gmsh.model.geo.add_line(p5, p2)

        # Right boundary
        self.c_right = gmsh.model.geo.add_line(p2, p3)

        # Top boundary (surface, in two parts)
        self.c_top2 = gmsh.model.geo.add_line(p3, p6)
        self.c_top1 = gmsh.model.geo.add_line(p6, p4)

        # Left boundary
        self.c_left = gmsh.model.geo.add_line(p4, p1)

        # Fault surface (from surface to fault tip)
        self.c_fault = gmsh.model.geo.add_line(p6, p7)

        # Extension below fault tip (not a fault, just mesh transition)
        self.c_fault_ext = gmsh.model.geo.add_line(p7, p5)

        # Create surfaces
        # Left block (negative x side)
        c_left_loop = gmsh.model.geo.add_curve_loop([
            self.c_bottom1, -self.c_fault_ext, -self.c_fault, self.c_top1, self.c_left
        ])
        self.s_left = gmsh.model.geo.add_plane_surface([c_left_loop])

        # Right block (positive x side)
        c_right_loop = gmsh.model.geo.add_curve_loop([
            self.c_bottom2, self.c_right, self.c_top2, self.c_fault, self.c_fault_ext
        ])
        self.s_right = gmsh.model.geo.add_plane_surface([c_right_loop])

        gmsh.model.geo.synchronize()

    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.
        """
        # Materials (both sides have same thermoelastic properties)
        materials = (
            MaterialGroup(tag=1, entities=[self.s_left]),
            MaterialGroup(tag=2, entities=[self.s_right]),
        )
        for material in materials:
            material.create_physical_group()

        # Boundary groups
        boundary_groups = (
            BoundaryGroup(name="boundary_left", tag=10, dim=1, entities=[self.c_left]),
            BoundaryGroup(name="boundary_right", tag=11, dim=1, entities=[self.c_right]),
            BoundaryGroup(name="boundary_bottom", tag=12, dim=1, entities=[self.c_bottom1, self.c_bottom2]),
            BoundaryGroup(name="boundary_top", tag=13, dim=1, entities=[self.c_top1, self.c_top2]),
            BoundaryGroup(name="fault", tag=20, dim=1, entities=[self.c_fault]),
        )
        for group in boundary_groups:
            group.create_physical_group()

        # Buried fault edge (vertex group for the fault tip)
        gmsh.model.geo.synchronize()
        # The fault tip is at the end of c_fault (point p7)
        fault_boundary = gmsh.model.getBoundary([(1, self.c_fault)], oriented=False)
        for dim, tag in fault_boundary:
            coord = gmsh.model.getValue(dim, tag, [])
            if coord[1] < -1.0:  # Below surface
                gmsh.model.addPhysicalGroup(0, [tag], 21)
                gmsh.model.setPhysicalName(0, 21, "fault_edge")
                break

    def generate_mesh(self, cell):
        """Generate the mesh.
        """
        # Disable default sizing
        gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)

        # Distance from fault for mesh sizing
        field_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(field_distance, "CurvesList", [self.c_fault])

        # Mesh size field with geometric progression from fault
        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(field_distance, min_dx=self.DX_FAULT, bias=self.DX_BIAS)
        gmsh.model.mesh.field.setString(field_size, "F", math_exp)

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


# End of file
