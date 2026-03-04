#!/usr/bin/env nemesis
"""Generate a tri mesh for the ALS model domain using Gmsh.

Fault curves are loaded from node coordinates extracted from the input mesh file
20260206_130505_pylith_model.msh. Two fault traces run diagonally (NW-SE) across
the rectangular domain.

Run `generate_gmsh.py --help` to see the command line options.
Run `generate_gmsh.py --write` to generate the mesh.
"""
import gmsh
from pylith.meshio.gmsh_utils import (VertexGroup, MaterialGroup, GenerateMesh)


class App(GenerateMesh):
    """
    Application for generating the mesh.

    Domain (meters):
        369921 <= x <= 388747
        6527123 <= y <= 6547901

    p4-----------------p3
    |                   |
    |  fault_a          |
    |     \    fault_b  |
    |      \    /       |
    |       \  /        |
    |        \/         |
    p1-----------------p2
    """
    X_WEST = 369921.0
    X_EAST = 388747.0
    Y_SOUTH = 6527123.0
    Y_NORTH = 6547901.0

    DX_FAULT = 200.0
    DX_BIAS = 1.05

    # Fault A: NW-SE trace (57 control points, extracted from input mesh)
    FAULT_A_POINTS = (
        (369980.4, 6542902.3),
        (370190.4, 6542655.3),
        (370352.4, 6542464.9),
        (370547.4, 6542210.7),
        (370699.6, 6542012.3),
        (370851.7, 6541813.9),
        (371003.9, 6541615.6),
        (371156.0, 6541417.2),
        (371416.5, 6541133.5),
        (371585.5, 6540949.3),
        (371754.6, 6540765.1),
        (371923.6, 6540580.9),
        (372092.7, 6540396.8),
        (372314.7, 6540079.2),
        (372458.0, 6539874.3),
        (372648.7, 6539595.9),
        (372790.0, 6539389.6),
        (373063.6, 6538971.2),
        (373200.4, 6538762.0),
        (373337.2, 6538552.7),
        (373474.0, 6538343.5),
        (373610.8, 6538134.2),
        (373737.0, 6537900.1),
        (373855.5, 6537680.0),
        (373974.1, 6537459.9),
        (374092.7, 6537239.8),
        (374211.3, 6537019.7),
        (374329.8, 6536799.6),
        (374494.1, 6536528.4),
        (374623.7, 6536314.6),
        (374753.2, 6536100.8),
        (374882.8, 6535887.0),
        (375012.3, 6535673.1),
        (375251.7, 6535539.7),
        (375470.1, 6535418.0),
        (375688.4, 6535296.3),
        (375906.8, 6535174.6),
        (376125.2, 6535052.8),
        (376343.6, 6534931.1),
        (376679.1, 6534813.6),
        (376956.1, 6534678.3),
        (377180.7, 6534568.5),
        (377405.3, 6534458.7),
        (377504.7, 6534410.1),
        (377670.8, 6534067.1),
        (377779.8, 6533842.1),
        (377888.8, 6533617.1),
        (377942.6, 6533561.1),
        (377962.7, 6533510.6),
        (378399.6, 6533388.1),
        (378797.6, 6533175.1),
        (379018.0, 6533057.1),
        (379167.2, 6533032.6),
        (379417.4, 6532884.0),
        (379632.3, 6532756.3),
        (379847.2, 6532628.6),
        (380062.1, 6532500.9),
    )

    # Fault B: N-SE trace (92 control points, extracted from input mesh)
    FAULT_B_POINTS = (
        (375550.4, 6547752.2),
        (375654.3, 6547496.2),
        (375785.6, 6547204.5),
        (375888.3, 6546976.6),
        (375990.9, 6546748.6),
        (376093.5, 6546520.6),
        (376196.1, 6546292.7),
        (376330.8, 6546012.1),
        (376438.9, 6545786.7),
        (376547.1, 6545561.3),
        (376655.3, 6545335.9),
        (376763.4, 6545110.5),
        (376871.6, 6544885.1),
        (376979.7, 6544659.7),
        (377126.6, 6544353.5),
        (377234.8, 6544128.1),
        (377342.9, 6543902.7),
        (377536.3, 6543540.2),
        (377654.0, 6543319.6),
        (377771.7, 6543099.1),
        (377889.4, 6542878.5),
        (378007.1, 6542657.9),
        (378124.8, 6542437.4),
        (378242.4, 6542216.8),
        (378360.1, 6541996.3),
        (378477.8, 6541775.7),
        (378710.2, 6541395.1),
        (378840.4, 6541181.8),
        (378970.7, 6540968.4),
        (379101.0, 6540755.0),
        (379231.2, 6540541.6),
        (379361.5, 6540328.3),
        (379491.8, 6540114.9),
        (379622.1, 6539901.5),
        (379752.3, 6539688.1),
        (379882.6, 6539474.8),
        (380012.9, 6539261.4),
        (380143.1, 6539048.0),
        (380273.4, 6538834.6),
        (380403.7, 6538621.3),
        (380534.0, 6538407.9),
        (380664.2, 6538194.5),
        (380794.5, 6537981.1),
        (380924.8, 6537767.8),
        (380975.5, 6537684.6),
        (381194.2, 6537385.4),
        (381341.7, 6537183.5),
        (381489.2, 6536981.7),
        (381636.6, 6536779.8),
        (381784.1, 6536578.0),
        (381931.6, 6536376.1),
        (382162.8, 6536095.2),
        (382321.7, 6535902.2),
        (382480.6, 6535709.2),
        (382648.8, 6535504.8),
        (382807.7, 6535311.8),
        (382966.5, 6535118.7),
        (383125.4, 6534925.7),
        (383284.2, 6534732.6),
        (383443.1, 6534539.6),
        (383625.0, 6534326.5),
        (383787.3, 6534136.4),
        (383949.7, 6533946.3),
        (384112.0, 6533756.1),
        (384274.3, 6533566.0),
        (384436.7, 6533375.9),
        (384599.0, 6533185.7),
        (384761.3, 6532995.6),
        (384923.6, 6532805.5),
        (385086.0, 6532615.3),
        (385248.3, 6532425.2),
        (385410.6, 6532235.1),
        (385628.8, 6532003.5),
        (385821.9, 6531798.5),
        (385993.3, 6531616.5),
        (386164.7, 6531434.6),
        (386411.5, 6531127.5),
        (386568.1, 6530932.6),
        (386724.8, 6530737.8),
        (386881.4, 6530542.9),
        (387038.0, 6530348.0),
        (387294.6, 6529980.1),
        (387437.6, 6529775.0),
        (387679.4, 6529428.0),
        (387885.3, 6529079.0),
        (388012.3, 6528863.7),
        (388139.3, 6528648.3),
        (388266.3, 6528433.0),
        (388370.3, 6528257.0),
        (388537.4, 6527892.1),
        (388641.4, 6527664.8),
        (388745.5, 6527437.5),
    )

    def __init__(self):
        self.cell_choices = {
            "default": "tri",
            "choices": ["tri"],
        }
        self.filename = "mesh_tri.msh"

    def create_geometry(self):
        """Create geometry.

        Rectangular domain with two embedded fault curves.
        """
        # Domain corners
        p1 = gmsh.model.geo.add_point(self.X_WEST, self.Y_SOUTH, 0.0)
        p2 = gmsh.model.geo.add_point(self.X_EAST, self.Y_SOUTH, 0.0)
        p3 = gmsh.model.geo.add_point(self.X_EAST, self.Y_NORTH, 0.0)
        p4 = gmsh.model.geo.add_point(self.X_WEST, self.Y_NORTH, 0.0)

        # Domain boundary curves
        self.c_south = gmsh.model.geo.add_line(p1, p2)
        self.c_east = gmsh.model.geo.add_line(p2, p3)
        self.c_north = gmsh.model.geo.add_line(p3, p4)
        self.c_west = gmsh.model.geo.add_line(p4, p1)

        # Create domain surface
        loop = gmsh.model.geo.add_curve_loop([self.c_south, self.c_east, self.c_north, self.c_west])
        self.s_domain = gmsh.model.geo.add_plane_surface([loop])

        # Create fault A spline from extracted points
        pts_a = []
        for x, y in self.FAULT_A_POINTS:
            pts_a.append(gmsh.model.geo.add_point(x, y, 0.0))
        self.c_fault_a = gmsh.model.geo.add_spline(pts_a)

        # Create fault B spline from extracted points
        pts_b = []
        for x, y in self.FAULT_B_POINTS:
            pts_b.append(gmsh.model.geo.add_point(x, y, 0.0))
        self.c_fault_b = gmsh.model.geo.add_spline(pts_b)

        gmsh.model.geo.synchronize()

        # Embed fault curves in the domain surface so the mesh conforms to them
        gmsh.model.mesh.embed(1, [self.c_fault_a, self.c_fault_b], 2, self.s_domain)

    def mark(self):
        """Mark geometry for materials, boundary conditions, faults."""
        # Single material for the whole domain
        materials = (
            MaterialGroup(tag=1, entities=[self.s_domain]),
        )
        for material in materials:
            material.create_physical_group()

        # Boundary and fault groups
        vertex_groups = (
            VertexGroup(name="boundary_south", tag=10, dim=1, entities=[self.c_south]),
            VertexGroup(name="boundary_east", tag=11, dim=1, entities=[self.c_east]),
            VertexGroup(name="boundary_north", tag=12, dim=1, entities=[self.c_north]),
            VertexGroup(name="boundary_west", tag=13, dim=1, entities=[self.c_west]),
            VertexGroup(name="fault_a", tag=20, dim=1, entities=[self.c_fault_a]),
            VertexGroup(name="fault_b", tag=21, dim=1, entities=[self.c_fault_b]),
        )
        for group in vertex_groups:
            group.create_physical_group()

    def generate_mesh(self, cell):
        """Generate the mesh."""
        # Disable default sizing
        gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)

        # Distance field from both faults
        field_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(field_distance, "CurvesList", [self.c_fault_a, self.c_fault_b])

        # Size field: refine near faults, coarsen away
        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(field_distance, min_dx=self.DX_FAULT, bias=self.DX_BIAS)
        gmsh.model.mesh.field.setString(field_size, "F", math_exp)

        gmsh.model.mesh.field.setAsBackgroundMesh(field_size)

        gmsh.model.mesh.generate(2)
        gmsh.model.mesh.optimize("Laplace2D")


if __name__ == "__main__":
    App().main()
