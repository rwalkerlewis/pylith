#!/usr/bin/env nemesis
"""Generate a 2D tri or quad mesh for the ALS paleoseismic example using Gmsh.

The mesh covers a rectangular domain based on the ALS geospatial data
(EPSG:28351) with a simplified diagonal fault trace splitting the domain
into two material regions.  We use the ``gmsh_utils`` module provided
with PyLith for helper functions.

Domain (in metres, EPSG 28351 – GDA94 / MGA zone 51):
    X : 370 000 – 389 000   (19 km)
    Y : 6 527 000 – 6 548 000   (21 km)

The fault trace is a simplified representation of the main NE–SW-trending
fault digitised from the ALS geospatial layers.  It enters through the
south boundary and exits through the north boundary, dividing the domain
into a western ("slab_west") and eastern ("slab_east") material block.

Run ``generate_gmsh.py --help`` to see the command-line options.
Run ``generate_gmsh.py --write`` to generate the mesh.
"""

import gmsh

# Try to import PyLith's gmsh_utils; fall back to a minimal shim so
# the script can also be run without a full PyLith installation.
try:
    from pylith.meshio.gmsh_utils import (
        BoundaryGroup,
        MaterialGroup,
        GenerateMesh,
    )
except ImportError:
    # Minimal fallback so that the mesh can still be generated
    # without a full PyLith install.
    import argparse, sys

    class MaterialGroup:
        def __init__(self, tag, entities):
            self.tag = tag
            self.entities = entities
        def create_physical_group(self):
            gmsh.model.addPhysicalGroup(2, self.entities, self.tag)
            gmsh.model.setPhysicalName(2, self.tag,
                                       f"material-id:{self.tag}")

    class BoundaryGroup:
        def __init__(self, name, tag, dim, entities):
            self.name = name
            self.tag = tag
            self.dim = dim
            self.entities = entities
        def create_physical_group(self):
            gmsh.model.addPhysicalGroup(self.dim, self.entities, self.tag)
            gmsh.model.setPhysicalName(self.dim, self.tag, self.name)

    class GenerateMesh:
        """Minimal base class when pylith is not available."""
        cell_choices = {"default": "tri", "choices": ["tri", "quad"]}
        filename = "mesh_tri.msh"

        @staticmethod
        def get_math_progression(field_distance, min_dx, bias):
            return f"{min_dx}*{bias}^(F{field_distance}/({min_dx}))"

        def create_geometry(self):
            raise NotImplementedError

        def mark(self):
            raise NotImplementedError

        def generate_mesh(self, cell):
            raise NotImplementedError

        def main(self):
            parser = argparse.ArgumentParser()
            parser.add_argument("--write", action="store_true")
            parser.add_argument("--cell", default="tri",
                                choices=["tri", "quad"])
            args = parser.parse_args()
            if not args.write:
                print("Use --write to generate mesh.")
                return
            gmsh.initialize()
            gmsh.model.add("als")
            self.create_geometry()
            self.mark()
            self.generate_mesh(args.cell)
            gmsh.write(self.filename)
            gmsh.finalize()


class App(GenerateMesh):
    """Application to generate the ALS paleoseismic mesh using Gmsh.

    Domain (m, EPSG 28351):
        X: 370 000 – 389 000
        Y: 6 527 000 – 6 548 000

    Simplified topology::

        p4----p_fn----p3
        |      /       |
        |     /        |
        |    /         |
        |   /  east    |
        |  /           |
        | / west       |
        p1----p_fs----p2

    The fault trace runs from p_fs on the south boundary to p_fn on the
    north boundary.
    """

    # ---- Domain extent (metres) -----------------------------------------------
    DOMAIN_W = 370000.0   # West boundary
    DOMAIN_E = 389000.0   # East boundary
    DOMAIN_S = 6527000.0  # South boundary
    DOMAIN_N = 6548000.0  # North boundary

    # ---- Fault trace key-points (simplified from geospatial data) -------------
    # The fault enters at the south boundary and exits at the north boundary.
    # These waypoints capture the main NE–SW trend observed in the data.
    FAULT_TRACE = [
        (388000.0, 6527000.0),   # south boundary intersection
        (387400.0, 6530000.0),
        (385500.0, 6532000.0),
        (383500.0, 6534500.0),
        (381900.0, 6536400.0),
        (380900.0, 6537700.0),
        (379100.0, 6540000.0),
        (377500.0, 6543500.0),
        (376200.0, 6546300.0),
        (375500.0, 6548000.0),   # north boundary intersection
    ]

    # ---- Mesh discretisation --------------------------------------------------
    DX = 1000.0        # nominal cell size (m)

    def __init__(self):
        self.cell_choices = {
            "default": "tri",
            "choices": ["tri", "quad"],
        }
        self.filename = "mesh_tri.msh"

    # ------------------------------------------------------------------ helpers
    @staticmethod
    def _add_points(coords):
        """Add Gmsh points and return their tags."""
        return [gmsh.model.geo.add_point(x, y, 0.0) for x, y in coords]

    # --------------------------------------------------------------- geometry
    def create_geometry(self):
        """Create geometry: rectangular domain split by a fault."""

        # Corner points (counter-clockwise from SW)
        p1 = gmsh.model.geo.add_point(self.DOMAIN_W, self.DOMAIN_S, 0.0)
        p2 = gmsh.model.geo.add_point(self.DOMAIN_E, self.DOMAIN_S, 0.0)
        p3 = gmsh.model.geo.add_point(self.DOMAIN_E, self.DOMAIN_N, 0.0)
        p4 = gmsh.model.geo.add_point(self.DOMAIN_W, self.DOMAIN_N, 0.0)

        # Fault end-points (on the boundary)
        ft = self.FAULT_TRACE
        p_fs = gmsh.model.geo.add_point(ft[0][0], ft[0][1], 0.0)   # south
        p_fn = gmsh.model.geo.add_point(ft[-1][0], ft[-1][1], 0.0)  # north

        # Interior fault waypoints
        fault_interior = self._add_points(ft[1:-1])
        fault_pts = [p_fs] + fault_interior + [p_fn]

        # ---- Boundary curves -------------------------------------------------
        # South boundary is split by the fault entry point:
        #   p1 --yneg_w--> p_fs --yneg_e--> p2
        self.c_yneg_w = gmsh.model.geo.add_line(p1, p_fs)
        self.c_yneg_e = gmsh.model.geo.add_line(p_fs, p2)

        # East boundary
        self.c_xpos = gmsh.model.geo.add_line(p2, p3)

        # North boundary is split by the fault exit point:
        #   p3 --ypos_e--> p_fn --ypos_w--> p4
        self.c_ypos_e = gmsh.model.geo.add_line(p3, p_fn)
        self.c_ypos_w = gmsh.model.geo.add_line(p_fn, p4)

        # West boundary
        self.c_xneg = gmsh.model.geo.add_line(p4, p1)

        # ---- Fault curve (polyline south → north) ----------------------------
        fault_segments = []
        for i in range(len(fault_pts) - 1):
            seg = gmsh.model.geo.add_line(fault_pts[i], fault_pts[i + 1])
            fault_segments.append(seg)
        self.c_fault_segments = fault_segments

        # ---- Surfaces ---------------------------------------------------------
        # West surface (left of fault, traversed CCW):
        #   yneg_w  →  fault  →  ypos_w  →  xneg
        west_loop = gmsh.model.geo.add_curve_loop(
            [self.c_yneg_w] + fault_segments + [self.c_ypos_w, self.c_xneg]
        )
        self.s_west = gmsh.model.geo.add_plane_surface([west_loop])

        # East surface (right of fault, traversed CCW):
        #   yneg_e  →  xpos  →  ypos_e  →  -fault
        east_loop = gmsh.model.geo.add_curve_loop(
            [self.c_yneg_e, self.c_xpos, self.c_ypos_e]
            + [-s for s in reversed(fault_segments)]
        )
        self.s_east = gmsh.model.geo.add_plane_surface([east_loop])

        gmsh.model.geo.synchronize()

    # ------------------------------------------------------------------- mark
    def mark(self):
        """Mark geometry for materials, boundaries, and faults."""

        # Material groups (tag must match the label_value in PyLith cfg)
        materials = (
            MaterialGroup(tag=1, entities=[self.s_west]),
            MaterialGroup(tag=2, entities=[self.s_east]),
        )
        for m in materials:
            m.create_physical_group()

        # Boundary groups
        boundary_groups = (
            BoundaryGroup(name="boundary_xneg", tag=10, dim=1,
                          entities=[self.c_xneg]),
            BoundaryGroup(name="boundary_xpos", tag=11, dim=1,
                          entities=[self.c_xpos]),
            BoundaryGroup(name="boundary_yneg", tag=12, dim=1,
                          entities=[self.c_yneg_w, self.c_yneg_e]),
            BoundaryGroup(name="boundary_ypos", tag=13, dim=1,
                          entities=[self.c_ypos_w, self.c_ypos_e]),
            BoundaryGroup(name="fault", tag=20, dim=1,
                          entities=self.c_fault_segments),
        )
        for bg in boundary_groups:
            bg.create_physical_group()

    # ------------------------------------------------------------- mesh gen
    def generate_mesh(self, cell):
        """Generate 2-D mesh.

        Args:
            cell: ``"tri"`` or ``"quad"``
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
