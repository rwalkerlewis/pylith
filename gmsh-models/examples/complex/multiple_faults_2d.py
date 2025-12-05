#!/usr/bin/env nemesis
"""Generate a 2D mesh with multiple intersecting faults using Gmsh.

This complex example demonstrates:
- Domain with multiple faults (main fault + two splay faults)
- Multiple materials created by fault intersections
- Multiple fault marking for PyLith
- Complex mesh refinement near multiple faults
- Fault intersection handling

Domain: 200km x 150km
-100.0 km <= x <= 100.0 km
-75.0 km <= y <= 75.0 km

Faults:
- Main fault: Vertical line at x=0
- West splay: Diagonal fault from main fault to west boundary
- East splay: Diagonal fault from main fault to east boundary

Run `multiple_faults_2d.py --write` to generate the mesh.
"""

import math
import gmsh
from pylith.meshio.gmsh_utils import (BoundaryGroup, MaterialGroup, GenerateMesh)


class App(GenerateMesh):
    """
    Application used to generate a 2D mesh with multiple intersecting faults.
    
    This example shows:
    - Multiple faults creating multiple material regions
    - Fault intersection points
    - Complex geometry handling
    - Refinement near multiple faults
    """
    
    DOMAIN_X = 200.0e+3  # 200 km
    DOMAIN_Y = 150.0e+3  # 150 km
    
    # Fault geometry parameters
    SPLAY_ANGLE = 30.0  # Angle of splay faults from vertical (degrees)
    SPLAY_Y_POS = 50.0e+3  # Y position where splays branch from main fault
    
    DX_FAULT = 3.0e+3    # Cell size on faults: 3 km
    DX_BIAS = 1.12       # Bias factor for cell size growth
    
    def __init__(self):
        """Constructor.
        """
        super().__init__()
        self.cell_choices = {
            "default": "tri",
            "choices": ["tri"],
        }
        self.filename = "multiple_faults_2d.msh"
    
    def create_geometry(self):
        """Create geometry with multiple intersecting faults.
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
        
        # Create points on main fault (vertical at x=0)
        p5 = gmsh.model.geo.add_point(0.0, y1, 0.0)
        p6 = gmsh.model.geo.add_point(0.0, y1 + ly, 0.0)
        self.p_intersect = gmsh.model.geo.add_point(0.0, self.SPLAY_Y_POS, 0.0)
        
        # Calculate splay fault endpoints
        splay_angle_rad = math.radians(self.SPLAY_ANGLE)
        splay_length = 0.5 * lx / math.cos(splay_angle_rad)
        
        # West splay endpoints
        p_west_bottom = gmsh.model.geo.add_point(
            -splay_length * math.sin(splay_angle_rad),
            y1,
            0.0
        )
        p_west_top = gmsh.model.geo.add_point(
            -splay_length * math.sin(splay_angle_rad),
            y1 + ly,
            0.0
        )
        
        # East splay endpoints
        p_east_bottom = gmsh.model.geo.add_point(
            splay_length * math.sin(splay_angle_rad),
            y1,
            0.0
        )
        p_east_top = gmsh.model.geo.add_point(
            splay_length * math.sin(splay_angle_rad),
            y1 + ly,
            0.0
        )
        
        # Create boundary curves
        self.c_yneg_west = gmsh.model.geo.add_line(p1, p_west_bottom)
        self.c_yneg_center = gmsh.model.geo.add_line(p_west_bottom, p_east_bottom)
        self.c_yneg_east = gmsh.model.geo.add_line(p_east_bottom, p2)
        self.c_xpos = gmsh.model.geo.add_line(p2, p3)
        self.c_ypos_east = gmsh.model.geo.add_line(p3, p_east_top)
        self.c_ypos_center = gmsh.model.geo.add_line(p_east_top, p_west_top)
        self.c_ypos_west = gmsh.model.geo.add_line(p_west_top, p4)
        self.c_xneg = gmsh.model.geo.add_line(p4, p1)
        
        # Create main fault segments
        self.c_fault_main_bottom = gmsh.model.geo.add_line(p5, self.p_intersect)
        self.c_fault_main_top = gmsh.model.geo.add_line(self.p_intersect, p6)
        
        # Create splay faults
        self.c_fault_west_bottom = gmsh.model.geo.add_line(p_west_bottom, self.p_intersect)
        self.c_fault_west_top = gmsh.model.geo.add_line(self.p_intersect, p_west_top)
        self.c_fault_east_bottom = gmsh.model.geo.add_line(p_east_bottom, self.p_intersect)
        self.c_fault_east_top = gmsh.model.geo.add_line(self.p_intersect, p_east_top)
        
        # Create surfaces (multiple regions created by faults)
        # Bottom-left region
        loop0 = gmsh.model.geo.add_curve_loop([
            self.c_yneg_west, -self.c_fault_west_bottom, 
            -self.c_fault_main_bottom, self.c_xneg
        ])
        self.s_bottom_left = gmsh.model.geo.add_plane_surface([loop0])
        
        # Bottom-center region
        loop1 = gmsh.model.geo.add_curve_loop([
            self.c_yneg_center, -self.c_fault_east_bottom,
            self.c_fault_main_bottom, self.c_fault_west_bottom
        ])
        self.s_bottom_center = gmsh.model.geo.add_plane_surface([loop1])
        
        # Bottom-right region
        loop2 = gmsh.model.geo.add_curve_loop([
            self.c_yneg_east, self.c_xpos, self.c_ypos_east,
            -self.c_fault_east_bottom
        ])
        self.s_bottom_right = gmsh.model.geo.add_plane_surface([loop2])
        
        # Top-left region
        loop3 = gmsh.model.geo.add_curve_loop([
            self.c_fault_main_top, self.c_fault_west_top,
            self.c_ypos_west, self.c_xneg
        ])
        self.s_top_left = gmsh.model.geo.add_plane_surface([loop3])
        
        # Top-center region
        loop4 = gmsh.model.geo.add_curve_loop([
            -self.c_fault_main_top, -self.c_fault_east_top,
            self.c_ypos_center, -self.c_fault_west_top
        ])
        self.s_top_center = gmsh.model.geo.add_plane_surface([loop4])
        
        # Top-right region
        loop5 = gmsh.model.geo.add_curve_loop([
            self.c_fault_east_top, self.c_ypos_east,
            -self.c_ypos_center
        ])
        self.s_top_right = gmsh.model.geo.add_plane_surface([loop5])
        
        gmsh.model.geo.synchronize()
    
    def mark(self):
        """Mark geometry for materials, boundary conditions, and faults.
        """
        # Create materials for each region
        materials = (
            MaterialGroup(tag=1, entities=[self.s_bottom_left]),
            MaterialGroup(tag=2, entities=[self.s_bottom_center]),
            MaterialGroup(tag=3, entities=[self.s_bottom_right]),
            MaterialGroup(tag=4, entities=[self.s_top_left]),
            MaterialGroup(tag=5, entities=[self.s_top_center]),
            MaterialGroup(tag=6, entities=[self.s_top_right]),
        )
        for material in materials:
            material.create_physical_group()
        
        # Create boundary groups
        boundary_groups = (
            BoundaryGroup(name="boundary_xneg", tag=10, dim=1, entities=[self.c_xneg]),
            BoundaryGroup(name="boundary_xpos", tag=11, dim=1, entities=[self.c_xpos]),
            BoundaryGroup(name="boundary_yneg", tag=12, dim=1, entities=[
                self.c_yneg_west, self.c_yneg_center, self.c_yneg_east
            ]),
            BoundaryGroup(name="boundary_ypos", tag=13, dim=1, entities=[
                self.c_ypos_west, self.c_ypos_center, self.c_ypos_east
            ]),
        )
        for group in boundary_groups:
            group.create_physical_group()
        
        # Create fault groups
        fault_groups = (
            BoundaryGroup(name="fault_main", tag=20, dim=1, entities=[
                self.c_fault_main_bottom, self.c_fault_main_top
            ]),
            BoundaryGroup(name="fault_west", tag=21, dim=1, entities=[
                self.c_fault_west_bottom, self.c_fault_west_top
            ]),
            BoundaryGroup(name="fault_east", tag=22, dim=1, entities=[
                self.c_fault_east_bottom, self.c_fault_east_top
            ]),
            BoundaryGroup(name="fault_intersection", tag=30, dim=0, entities=[self.p_intersect]),
        )
        for group in fault_groups:
            group.create_physical_group()
    
    def generate_mesh(self, cell):
        """Generate mesh with refinement near all faults.
        """
        # Turn off default sizing methods
        gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)
        
        # Create distance field from all faults
        all_fault_curves = [
            self.c_fault_main_bottom, self.c_fault_main_top,
            self.c_fault_west_bottom, self.c_fault_west_top,
            self.c_fault_east_bottom, self.c_fault_east_top
        ]
        field_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(field_distance, "CurvesList", all_fault_curves)
        
        # Create size field with geometric progression
        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(field_distance, min_dx=self.DX_FAULT, bias=self.DX_BIAS)
        gmsh.model.mesh.field.setString(field_size, "F", math_exp)
        
        # Use the size field as background mesh
        gmsh.model.mesh.field.setAsBackgroundMesh(field_size)
        
        gmsh.model.mesh.generate(2)
        gmsh.model.mesh.optimize("Laplace2D")


if __name__ == "__main__":
    App().main()
