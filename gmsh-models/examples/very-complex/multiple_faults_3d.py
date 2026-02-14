#!/usr/bin/env nemesis
"""Generate a 3D mesh with multiple intersecting faults using Gmsh.

This very complex example demonstrates:
- 3D domain with multiple faults
- Multiple materials created by fault intersections in 3D
- Multiple fault marking for PyLith in 3D
- Complex 3D mesh refinement near multiple faults
- 3D fault intersection handling

Domain: 100km x 100km x 50km (depth)
-50.0 km <= x <= 50.0 km
-50.0 km <= y <= 50.0 km
-50.0 km <= z <= 0.0 km (z=0 is surface)

Faults:
- Main fault: Vertical plane at x=0 (strike-slip)
- Secondary fault: Dipping plane intersecting main fault

Run `multiple_faults_3d.py --write` to generate the mesh.
"""

import math
import gmsh
from pylith.meshio.gmsh_utils import (BoundaryGroup, MaterialGroup, GenerateMesh)


class App(GenerateMesh):
    """
    Application used to generate a 3D mesh with multiple intersecting faults.
    
    This example shows:
    - 3D domain with multiple faults
    - Multiple materials in 3D
    - 3D fault marking
    - Complex 3D mesh refinement
    """
    
    DOMAIN_X = 100.0e+3  # 100 km
    DOMAIN_Y = 100.0e+3  # 100 km
    DOMAIN_Z = 50.0e+3   # 50 km depth
    
    # Fault geometry parameters
    FAULT_DIP = 45.0  # Dip angle of secondary fault (degrees)
    FAULT_STRIKE_OFFSET = 20.0e+3  # Offset of secondary fault along y-axis
    
    DX_FAULT = 2.0e+3    # Cell size on faults: 2 km
    DX_BIAS = 1.1        # Bias factor for cell size growth
    
    def __init__(self):
        """Constructor.
        """
        super().__init__()
        self.cell_choices = {
            "default": "tet",
            "choices": ["tet"],
        }
        self.filename = "multiple_faults_3d.msh"
    
    def create_geometry(self):
        """Create 3D geometry with multiple intersecting faults using OpenCASCADE.
        """
        lx = self.DOMAIN_X
        ly = self.DOMAIN_Y
        lz = self.DOMAIN_Z
        x0 = -0.5 * lx
        y0 = -0.5 * ly
        z0 = -lz  # Bottom of domain
        
        # Create main domain box
        self.v_domain = gmsh.model.occ.add_box(x0, y0, z0, lx, ly, lz)
        
        # Create main fault plane (vertical plane at x=0)
        # Create a large rectangle that will be used to cut the domain
        fault_thickness = 1.0e+3  # Small thickness for the fault plane
        self.v_fault_main = gmsh.model.occ.add_box(
            -fault_thickness/2, y0, z0,
            fault_thickness, ly, lz
        )
        
        # Create secondary fault plane (dipping plane)
        fault_dip_rad = math.radians(self.FAULT_DIP)
        y_offset = self.FAULT_STRIKE_OFFSET
        
        # Create a box for the secondary fault, then rotate it
        fault_width = ly * 1.5  # Make sure it spans the domain
        fault_height = lz / math.cos(fault_dip_rad)  # Account for dip
        self.v_fault_sec_temp = gmsh.model.occ.add_box(
            0.0, y_offset, z0,
            lx, fault_thickness, fault_height
        )
        
        # Rotate secondary fault to create dip
        gmsh.model.occ.rotate(
            [(3, self.v_fault_sec_temp)],
            0.0, y_offset, z0,
            0.0, 1.0, 0.0,
            fault_dip_rad
        )
        self.v_fault_sec = self.v_fault_sec_temp
        
        # Synchronize geometry
        gmsh.model.occ.synchronize()
        
        # Fragment the domain with the faults to create separate volumes
        all_objects, all_maps = gmsh.model.occ.fragment(
            [(3, self.v_domain)],
            [(3, self.v_fault_main), (3, self.v_fault_sec)]
        )
        
        gmsh.model.occ.synchronize()
        
        # Extract volumes (domains) and surfaces (faults)
        self.volumes = [v[1] for v in all_objects if v[0] == 3]
        
        # Get fault surfaces from the mapping
        # The fault volumes become surfaces after fragmenting
        fault_surfaces = []
        for dim, tag in all_maps[1]:  # Maps from fault volumes
            if dim == 2:  # Surfaces
                fault_surfaces.append(tag)
        
        # Identify main and secondary fault surfaces by position
        # Main fault should be near x=0, secondary fault should be at y_offset
        self.s_fault_main = None
        self.s_fault_sec = None
        
        for surf_tag in fault_surfaces:
            bbox = gmsh.model.get_bounding_box(2, surf_tag)
            center_x = (bbox[0] + bbox[3]) / 2.0
            center_y = (bbox[1] + bbox[4]) / 2.0
            
            if abs(center_x) < lx * 0.1:  # Near x=0
                self.s_fault_main = surf_tag
            elif abs(center_y - y_offset) < ly * 0.1:  # Near y_offset
                self.s_fault_sec = surf_tag
        
        # Get boundary surfaces using bounding boxes
        dx = 100.0
        bbox = gmsh.model.get_bounding_box(-1, -1)
        get_bbox_entities = gmsh.model.get_entities_in_bounding_box
        
        self.s_xneg = get_bbox_entities(
            bbox[0]-dx, bbox[1]-dx, bbox[2]-dx,
            bbox[0]+dx, bbox[4]+dx, bbox[5]+dx, dim=2
        )
        self.s_xpos = get_bbox_entities(
            bbox[3]-dx, bbox[1]-dx, bbox[2]-dx,
            bbox[3]+dx, bbox[4]+dx, bbox[5]+dx, dim=2
        )
        self.s_yneg = get_bbox_entities(
            bbox[0]-dx, bbox[1]-dx, bbox[2]-dx,
            bbox[3]+dx, bbox[1]+dx, bbox[5]+dx, dim=2
        )
        self.s_ypos = get_bbox_entities(
            bbox[0]-dx, bbox[4]-dx, bbox[2]-dx,
            bbox[3]+dx, bbox[4]+dx, bbox[5]+dx, dim=2
        )
        self.s_zneg = get_bbox_entities(
            bbox[0]-dx, bbox[1]-dx, bbox[2]-dx,
            bbox[3]+dx, bbox[4]+dx, bbox[2]+dx, dim=2
        )
        self.s_zpos = get_bbox_entities(
            bbox[0]-dx, bbox[1]-dx, bbox[5]-dx,
            bbox[3]+dx, bbox[4]+dx, bbox[5]+dx, dim=2
        )
    
    def mark(self):
        """Mark geometry for materials, boundary conditions, and faults.
        """
        # Create materials for each volume
        # In a real scenario, you'd identify which volumes correspond to which materials
        # For simplicity, we'll assign materials sequentially
        materials = []
        for i, vol in enumerate(self.volumes, start=1):
            materials.append(MaterialGroup(tag=i, entities=[vol]))
        
        for material in materials:
            material.create_physical_group()
        
        # Create boundary groups
        boundary_groups = (
            BoundaryGroup(name="boundary_xneg", tag=10, dim=2, entities=[tag for dim, tag in self.s_xneg]),
            BoundaryGroup(name="boundary_xpos", tag=11, dim=2, entities=[tag for dim, tag in self.s_xpos]),
            BoundaryGroup(name="boundary_yneg", tag=12, dim=2, entities=[tag for dim, tag in self.s_yneg]),
            BoundaryGroup(name="boundary_ypos", tag=13, dim=2, entities=[tag for dim, tag in self.s_ypos]),
            BoundaryGroup(name="boundary_bottom", tag=14, dim=2, entities=[tag for dim, tag in self.s_zneg]),
            BoundaryGroup(name="boundary_top", tag=15, dim=2, entities=[tag for dim, tag in self.s_zpos]),
        )
        for group in boundary_groups:
            group.create_physical_group()
        
        # Create fault groups
        fault_surfaces = []
        if self.s_fault_main is not None:
            fault_surfaces.append(self.s_fault_main)
        if self.s_fault_sec is not None:
            fault_surfaces.append(self.s_fault_sec)
        
        if len(fault_surfaces) == 2:
            fault_groups = (
                BoundaryGroup(name="fault_main", tag=20, dim=2, entities=[self.s_fault_main]),
                BoundaryGroup(name="fault_secondary", tag=21, dim=2, entities=[self.s_fault_sec]),
            )
            for group in fault_groups:
                group.create_physical_group()
        elif len(fault_surfaces) == 1:
            # Fallback if only one fault was identified
            BoundaryGroup(name="fault", tag=20, dim=2, entities=fault_surfaces).create_physical_group()
    
    def generate_mesh(self, cell):
        """Generate 3D mesh with refinement near faults.
        """
        # Turn off default sizing methods
        gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)
        
        # Create distance field from all faults
        fault_surfaces = []
        if self.s_fault_main is not None:
            fault_surfaces.append(self.s_fault_main)
        if self.s_fault_sec is not None:
            fault_surfaces.append(self.s_fault_sec)
        
        if fault_surfaces:
            field_distance = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(field_distance, "SurfacesList", fault_surfaces)
            
            # Create size field with geometric progression
            field_size = gmsh.model.mesh.field.add("MathEval")
            math_exp = GenerateMesh.get_math_progression(field_distance, min_dx=self.DX_FAULT, bias=self.DX_BIAS)
            gmsh.model.mesh.field.setString(field_size, "F", math_exp)
            
            # Use the size field as background mesh
            gmsh.model.mesh.field.setAsBackgroundMesh(field_size)
        else:
            # Fallback to uniform mesh size if no faults identified
            gmsh.option.set_number("Mesh.MeshSizeMax", self.DX_FAULT * 5.0)
            gmsh.option.set_number("Mesh.MeshSizeMin", self.DX_FAULT * 5.0)
        
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Netgen")


if __name__ == "__main__":
    App().main()
