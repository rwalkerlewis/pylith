#!/usr/bin/env python3
"""
Generate a PyLith-compatible binary Gmsh mesh from point cloud data.

This script reads the original 20260206_190857_pylith_model.msh file
which contains 6815 points with material and fault assignments, performs
Delaunay triangulation, and writes a binary Gmsh 4.1 format mesh suitable
for PyLith simulations.

Materials (physical tags):
  2: unassigned
  3: A_b_YEG
  4: A_f_YEG
  5: A_me_st
  6: A_mgss_Y
  7: A_o_YEG
  8: A_s_YEG
  9: A_u_YEG

Fault (physical tag 1): fault_1
"""

import argparse
import struct
import numpy as np
from scipy.spatial import Delaunay
from pathlib import Path
from collections import Counter


# Physical group definitions matching the original mesh
PHYSICAL_GROUPS = {
    1: ("fault_1", 1),      # (name, dimension) - fault is 1D (edges)
    2: ("material-id:2", 2),  # 2D materials (triangles)
    3: ("material-id:3", 2),
    4: ("material-id:4", 2),
    5: ("material-id:5", 2),
    6: ("material-id:6", 2),
    7: ("material-id:7", 2),
    8: ("material-id:8", 2),
    9: ("material-id:9", 2),
}


def parse_msh_file(filepath):
    """
    Parse Gmsh 4.1 ASCII format file to extract points and their physical tags.
    
    Returns:
        points: numpy array of shape (N, 2) with x, y coordinates
        point_materials: numpy array of shape (N,) with material tag for each point
        fault_nodes: set of point indices that are on the fault
    """
    points = []
    point_materials = []
    fault_nodes = set()
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    lines = content.split('\n')
    i = 0
    
    while i < len(lines):
        line = lines[i].strip()
        
        if line == '$Entities':
            i += 1
            # Header: numPoints numCurves numSurfaces numVolumes
            header = lines[i].strip().split()
            num_points = int(header[0])
            i += 1
            
            for _ in range(num_points):
                parts = lines[i].strip().split()
                if len(parts) >= 5:
                    # Format: pointTag x y z numPhysicalTags [physicalTags...]
                    x = float(parts[1])
                    y = float(parts[2])
                    
                    num_phys = int(parts[4])
                    phys_tags = [int(parts[5 + j]) for j in range(num_phys)]
                    
                    point_idx = len(points)
                    points.append([x, y])
                    
                    # Determine material (highest non-fault tag, or 2 for unassigned)
                    material = 2  # default unassigned
                    for tag in phys_tags:
                        if tag == 1:
                            fault_nodes.add(point_idx)
                        elif 2 <= tag <= 9:
                            material = tag
                    
                    point_materials.append(material)
                i += 1
            
            # Skip curves, surfaces, volumes (not present in point-only mesh)
            while i < len(lines) and lines[i].strip() != '$EndEntities':
                i += 1
        
        i += 1
    
    return np.array(points), np.array(point_materials), fault_nodes


def triangulate_points(points):
    """
    Perform Delaunay triangulation on 2D points.
    
    Returns:
        triangles: numpy array of shape (M, 3) with vertex indices
    """
    tri = Delaunay(points)
    return tri.simplices


def compute_triangle_materials(triangles, point_materials):
    """
    Assign material to each triangle based on vertex materials.
    
    Uses majority vote among vertices. If tie, uses lowest material tag.
    
    Returns:
        triangle_materials: numpy array of shape (M,) with material tag
    """
    n_triangles = len(triangles)
    triangle_materials = np.zeros(n_triangles, dtype=int)
    
    for i, tri in enumerate(triangles):
        vertex_mats = point_materials[tri]
        counter = Counter(vertex_mats)
        # Get most common material (ties broken by lowest tag)
        most_common = counter.most_common()
        max_count = most_common[0][1]
        candidates = [mat for mat, count in most_common if count == max_count]
        triangle_materials[i] = min(candidates)
    
    return triangle_materials


def find_fault_edges(triangles, fault_nodes):
    """
    Find edges where both endpoints are fault nodes.
    
    Returns:
        fault_edges: list of (node1, node2) tuples
    """
    fault_edges = set()
    
    for tri in triangles:
        edges = [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])]
        for n1, n2 in edges:
            if n1 in fault_nodes and n2 in fault_nodes:
                # Ensure consistent ordering
                edge = (min(n1, n2), max(n1, n2))
                fault_edges.add(edge)
    
    return list(fault_edges)


def write_binary_msh(filepath, points, triangles, triangle_materials, 
                     fault_edges, physical_groups):
    """
    Write mesh in Gmsh 4.1 binary format.
    
    Binary format specs from Gmsh documentation:
    - size_t is 8 bytes (specified in header)
    - int is 4 bytes
    - double is 8 bytes
    """
    print(f"Writing binary mesh to {filepath}")
    
    n_points = len(points)
    n_triangles = len(triangles)
    n_fault_edges = len(fault_edges)
    
    # Collect which physical groups are actually used
    used_groups = set()
    if fault_edges:
        used_groups.add(1)  # fault
    used_groups.update(np.unique(triangle_materials))
    
    # Filter physical groups
    active_groups = {k: v for k, v in physical_groups.items() if k in used_groups}
    
    with open(filepath, 'wb') as f:
        # $MeshFormat - header is ASCII, then one binary int for endianness
        f.write(b"$MeshFormat\n")
        f.write(b"4.1 1 8\n")  # version, binary=1, sizeof(size_t)=8
        f.write(struct.pack('<i', 1))  # endian check (binary)
        f.write(b"\n$EndMeshFormat\n")
        
        # $PhysicalNames - entirely ASCII section
        f.write(b"$PhysicalNames\n")
        f.write(f"{len(active_groups)}\n".encode())
        for tag in sorted(active_groups.keys()):
            name, dim = active_groups[tag]
            f.write(f'{dim} {tag} "{name}"\n'.encode())
        f.write(b"$EndPhysicalNames\n")
        
        # $Entities - ASCII header line, then binary entity data
        n_curves = 1 if fault_edges else 0
        n_surfaces = len([t for t in active_groups if active_groups[t][1] == 2])
        
        f.write(b"$Entities\n")
        # numPoints numCurves numSurfaces numVolumes (size_t each, binary)
        f.write(struct.pack('<4Q', 0, n_curves, n_surfaces, 0))
        
        # Curve entity for fault (binary)
        # Format: curveTag(int) minX minY minZ maxX maxY maxZ(double*6) 
        #         numPhysicalTags(size_t) physicalTags[](int*) 
        #         numBoundingPoints(size_t) pointTags[](int*)
        if fault_edges:
            fault_pts = points[list(set(sum(fault_edges, ())))]
            min_x, min_y = fault_pts.min(axis=0)
            max_x, max_y = fault_pts.max(axis=0)
            f.write(struct.pack('<i', 1))  # curveTag
            f.write(struct.pack('<6d', min_x, min_y, 0.0, max_x, max_y, 0.0))  # bbox
            f.write(struct.pack('<Q', 1))  # numPhysicalTags (size_t)
            f.write(struct.pack('<i', 1))  # physicalTag (int)
            f.write(struct.pack('<Q', 0))  # numBoundingPoints (size_t)
        
        # Surface entities for materials (binary)
        # Format: surfaceTag(int) minX minY minZ maxX maxY maxZ(double*6)
        #         numPhysicalTags(size_t) physicalTags[](int*)
        #         numBoundingCurves(size_t) curveTags[](int*)
        for mat_tag in sorted([t for t in active_groups if active_groups[t][1] == 2]):
            mat_tris = triangles[triangle_materials == mat_tag]
            if len(mat_tris) > 0:
                mat_pts = points[np.unique(mat_tris.flatten())]
                min_x, min_y = mat_pts.min(axis=0)
                max_x, max_y = mat_pts.max(axis=0)
            else:
                min_x = min_y = max_x = max_y = 0.0
            f.write(struct.pack('<i', int(mat_tag)))  # surfaceTag (int)
            f.write(struct.pack('<6d', min_x, min_y, 0.0, max_x, max_y, 0.0))  # bbox
            f.write(struct.pack('<Q', 1))  # numPhysicalTags (size_t)
            f.write(struct.pack('<i', int(mat_tag)))  # physicalTag (int)
            f.write(struct.pack('<Q', 0))  # numBoundingCurves (size_t)
        
        f.write(b"\n$EndEntities\n")
        
        # $Nodes - binary
        # Format: numEntityBlocks(size_t) numNodes(size_t) minNodeTag(size_t) maxNodeTag(size_t)
        # Then for each entity block:
        #   entityDim(int) entityTag(int) parametric(int) numNodesInBlock(size_t)
        #   nodeTags[](size_t) coordinates[](double*3)
        f.write(b"$Nodes\n")
        num_entity_blocks = 1
        f.write(struct.pack('<4Q', num_entity_blocks, n_points, 1, n_points))
        
        # Single entity block containing all nodes
        f.write(struct.pack('<3i', 2, 1, 0))  # entityDim=2, entityTag=1, parametric=0
        f.write(struct.pack('<Q', n_points))  # numNodesInBlock (size_t)
        
        # Node tags (size_t each)
        for i in range(1, n_points + 1):
            f.write(struct.pack('<Q', i))
        
        # Node coordinates (3 doubles each)
        for pt in points:
            f.write(struct.pack('<3d', pt[0], pt[1], 0.0))
        
        f.write(b"\n$EndNodes\n")
        
        # $Elements - binary
        # Format: numEntityBlocks(size_t) numElements(size_t) minElementTag(size_t) maxElementTag(size_t)
        # Then for each entity block:
        #   entityDim(int) entityTag(int) elementType(int) numElementsInBlock(size_t)
        #   elementTag(size_t) nodeTags[](size_t) for each element
        f.write(b"$Elements\n")
        
        tri_by_mat = {}
        for i, mat in enumerate(triangle_materials):
            mat_int = int(mat)
            if mat_int not in tri_by_mat:
                tri_by_mat[mat_int] = []
            tri_by_mat[mat_int].append(i)
        
        num_entity_blocks = len(tri_by_mat) + (1 if fault_edges else 0)
        total_elements = n_triangles + n_fault_edges
        
        f.write(struct.pack('<4Q', num_entity_blocks, total_elements, 1, total_elements))
        
        element_tag = 1
        
        # Fault edges (1D elements, element type 1 = 2-node line)
        if fault_edges:
            f.write(struct.pack('<3i', 1, 1, 1))  # entityDim=1, entityTag=1, elementType=1
            f.write(struct.pack('<Q', n_fault_edges))  # numElementsInBlock
            
            for n1, n2 in fault_edges:
                f.write(struct.pack('<Q', element_tag))  # elementTag
                f.write(struct.pack('<2Q', n1 + 1, n2 + 1))  # nodeTags (1-indexed)
                element_tag += 1
        
        # Triangles (2D elements, element type 2 = 3-node triangle)
        for mat_tag in sorted(tri_by_mat.keys()):
            tri_indices = tri_by_mat[mat_tag]
            f.write(struct.pack('<3i', 2, mat_tag, 2))  # entityDim=2, entityTag, elementType=2
            f.write(struct.pack('<Q', len(tri_indices)))  # numElementsInBlock
            
            for idx in tri_indices:
                tri = triangles[idx]
                f.write(struct.pack('<Q', element_tag))  # elementTag
                f.write(struct.pack('<3Q', tri[0] + 1, tri[1] + 1, tri[2] + 1))  # nodeTags
                element_tag += 1
        
        f.write(b"\n$EndElements\n")
    
    print(f"  Nodes: {n_points}")
    print(f"  Triangles: {n_triangles}")
    print(f"  Fault edges: {n_fault_edges}")
    print(f"  Materials: {sorted(tri_by_mat.keys())}")


def write_ascii_msh(filepath, points, triangles, triangle_materials,
                    fault_edges, physical_groups):
    """
    Write mesh in Gmsh 4.1 ASCII format for debugging.
    """
    print(f"Writing ASCII mesh to {filepath}")
    
    n_points = len(points)
    n_triangles = len(triangles)
    n_fault_edges = len(fault_edges)
    
    # Collect which physical groups are actually used
    used_groups = set()
    if fault_edges:
        used_groups.add(1)
    used_groups.update(np.unique(triangle_materials))
    
    active_groups = {k: v for k, v in physical_groups.items() if k in used_groups}
    
    with open(filepath, 'w') as f:
        # MeshFormat
        f.write("$MeshFormat\n")
        f.write("4.1 0 8\n")  # version, ascii, sizeof(size_t)
        f.write("$EndMeshFormat\n")
        
        # PhysicalNames
        f.write("$PhysicalNames\n")
        f.write(f"{len(active_groups)}\n")
        for tag in sorted(active_groups.keys()):
            name, dim = active_groups[tag]
            f.write(f'{dim} {tag} "{name}"\n')
        f.write("$EndPhysicalNames\n")
        
        # Entities
        n_curves = 1 if fault_edges else 0
        n_surfaces = len([t for t in active_groups if active_groups[t][1] == 2])
        
        f.write("$Entities\n")
        f.write(f"0 {n_curves} {n_surfaces} 0\n")
        
        if fault_edges:
            fault_pts = points[list(set(sum(fault_edges, ())))]
            min_x, min_y = fault_pts.min(axis=0)
            max_x, max_y = fault_pts.max(axis=0)
            f.write(f"1 {min_x} {min_y} 0 {max_x} {max_y} 0 1 1 0\n")
        
        for mat_tag in sorted([t for t in active_groups if active_groups[t][1] == 2]):
            mat_tris = triangles[triangle_materials == mat_tag]
            if len(mat_tris) > 0:
                mat_pts = points[np.unique(mat_tris.flatten())]
                min_x, min_y = mat_pts.min(axis=0)
                max_x, max_y = mat_pts.max(axis=0)
            else:
                min_x = min_y = max_x = max_y = 0.0
            f.write(f"{mat_tag} {min_x} {min_y} 0 {max_x} {max_y} 0 1 {mat_tag} 0\n")
        
        f.write("$EndEntities\n")
        
        # Nodes
        f.write("$Nodes\n")
        f.write(f"1 {n_points} 1 {n_points}\n")
        f.write(f"2 1 0 {n_points}\n")  # dim=2, entityTag=1, parametric=0
        for i in range(1, n_points + 1):
            f.write(f"{i}\n")
        for pt in points:
            f.write(f"{pt[0]:.10f} {pt[1]:.10f} 0.0\n")
        f.write("$EndNodes\n")
        
        # Elements
        tri_by_mat = {}
        for i, mat in enumerate(triangle_materials):
            if mat not in tri_by_mat:
                tri_by_mat[mat] = []
            tri_by_mat[mat].append(i)
        
        num_entity_blocks = len(tri_by_mat) + (1 if fault_edges else 0)
        total_elements = n_triangles + n_fault_edges
        
        f.write("$Elements\n")
        f.write(f"{num_entity_blocks} {total_elements} 1 {total_elements}\n")
        
        element_tag = 1
        
        if fault_edges:
            f.write(f"1 1 1 {n_fault_edges}\n")  # dim=1, entityTag=1, type=1
            for n1, n2 in fault_edges:
                f.write(f"{element_tag} {n1 + 1} {n2 + 1}\n")
                element_tag += 1
        
        for mat_tag in sorted(tri_by_mat.keys()):
            tri_indices = tri_by_mat[mat_tag]
            f.write(f"2 {mat_tag} 2 {len(tri_indices)}\n")  # dim=2, entityTag, type=2
            for idx in tri_indices:
                tri = triangles[idx]
                f.write(f"{element_tag} {tri[0] + 1} {tri[1] + 1} {tri[2] + 1}\n")
                element_tag += 1
        
        f.write("$EndElements\n")
    
    print(f"  Nodes: {n_points}")
    print(f"  Triangles: {n_triangles}")
    print(f"  Fault edges: {n_fault_edges}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate PyLith-compatible mesh from point cloud data"
    )
    parser.add_argument(
        "--input", "-i",
        default="20260206_190857_pylith_model.msh",
        help="Input Gmsh file with point entities (default: 20260206_190857_pylith_model.msh)"
    )
    parser.add_argument(
        "--output", "-o",
        default="pylith_mesh.msh",
        help="Output mesh file (default: pylith_mesh.msh)"
    )
    parser.add_argument(
        "--ascii", "-a",
        action="store_true",
        help="Write ASCII format instead of binary"
    )
    parser.add_argument(
        "--debug", "-d",
        action="store_true",
        help="Write both ASCII and binary for debugging"
    )
    
    args = parser.parse_args()
    
    input_path = Path(args.input)
    if not input_path.is_absolute():
        input_path = Path(__file__).parent / input_path
    
    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = Path(__file__).parent / output_path
    
    print(f"Reading input mesh: {input_path}")
    points, point_materials, fault_nodes = parse_msh_file(input_path)
    
    print(f"  Points: {len(points)}")
    print(f"  Fault nodes: {len(fault_nodes)}")
    for mat in range(2, 10):
        count = np.sum(point_materials == mat)
        if count > 0:
            print(f"  Material {mat}: {count} points")
    
    print("Triangulating...")
    triangles = triangulate_points(points)
    print(f"  Created {len(triangles)} triangles")
    
    print("Assigning materials to triangles...")
    triangle_materials = compute_triangle_materials(triangles, point_materials)
    
    print("Finding fault edges...")
    fault_edges = find_fault_edges(triangles, fault_nodes)
    print(f"  Found {len(fault_edges)} fault edges")
    
    if args.ascii:
        write_ascii_msh(output_path, points, triangles, triangle_materials,
                       fault_edges, PHYSICAL_GROUPS)
    elif args.debug:
        # Write both formats
        ascii_path = output_path.with_suffix('.ascii.msh')
        write_ascii_msh(ascii_path, points, triangles, triangle_materials,
                       fault_edges, PHYSICAL_GROUPS)
        write_binary_msh(output_path, points, triangles, triangle_materials,
                        fault_edges, PHYSICAL_GROUPS)
    else:
        write_binary_msh(output_path, points, triangles, triangle_materials,
                        fault_edges, PHYSICAL_GROUPS)
    
    print("Done!")


if __name__ == "__main__":
    main()
