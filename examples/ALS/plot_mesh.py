#!/usr/bin/env python3
"""
Plot mesh elements colored by material section.

Usage:
    python plot_mesh.py mesh_tri.msh
    python plot_mesh.py mesh_quad.msh
    python plot_mesh.py --help
"""

import argparse
import gmsh
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch
import matplotlib.colors as mcolors


def get_material_colors(n_materials):
    """Generate distinct colors for materials."""
    if n_materials <= 10:
        colors = plt.cm.tab10.colors[:n_materials]
    else:
        colors = plt.cm.tab20.colors[:n_materials]
    return list(colors)


def plot_mesh(mesh_file, output_file=None, show_boundaries=True, show_fault=True):
    """
    Plot mesh elements colored by material section.
    
    Args:
        mesh_file: Path to the .msh file
        output_file: Output image file (default: mesh_file with .png extension)
        show_boundaries: Whether to highlight boundary edges
        show_fault: Whether to highlight fault edges
    """
    if output_file is None:
        output_file = mesh_file.rsplit('.', 1)[0] + '_materials.png'
    
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)  # Suppress output
    gmsh.open(mesh_file)
    
    # Get all nodes
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    coords = node_coords.reshape(-1, 3)
    node_map = {tag: i for i, tag in enumerate(node_tags)}
    
    # Get physical groups
    physical_groups_2d = gmsh.model.getPhysicalGroups(dim=2)
    physical_groups_1d = gmsh.model.getPhysicalGroups(dim=1)
    
    # Identify materials (2D physical groups)
    materials = []
    for dim, tag in physical_groups_2d:
        name = gmsh.model.getPhysicalName(dim, tag)
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        materials.append({
            'tag': tag,
            'name': name,
            'entities': entities
        })
    
    # Generate colors for materials
    material_colors = get_material_colors(len(materials))
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 11))
    
    legend_patches = []
    
    # Plot elements for each material
    for mat_idx, material in enumerate(materials):
        color = material_colors[mat_idx]
        mat_name = material['name'] if material['name'] else f"Material {material['tag']}"
        
        all_polygons = []
        
        for entity in material['entities']:
            elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(dim=2, tag=entity)
            
            for etype, etags, enodes in zip(elem_types, elem_tags, elem_node_tags):
                if etype == 2:  # Triangle
                    nodes_per_elem = 3
                elif etype == 3:  # Quad
                    nodes_per_elem = 4
                else:
                    continue
                
                n_elems = len(etags)
                enodes = enodes.reshape(n_elems, nodes_per_elem)
                
                for elem_nodes in enodes:
                    verts = [coords[node_map[n]][:2] for n in elem_nodes]
                    all_polygons.append(verts)
        
        if all_polygons:
            pc = PolyCollection(
                all_polygons, 
                facecolor=color, 
                edgecolor='black',
                linewidth=0.2, 
                alpha=0.8
            )
            ax.add_collection(pc)
            legend_patches.append(Patch(facecolor=color, edgecolor='black', 
                                       label=f'{mat_name} ({len(all_polygons)} elements)'))
    
    # Plot boundaries and faults
    boundary_lines = []
    fault_lines = []
    
    for dim, tag in physical_groups_1d:
        name = gmsh.model.getPhysicalName(dim, tag)
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        
        for entity in entities:
            etypes, etags, enodes = gmsh.model.mesh.getElements(dim=1, tag=entity)
            if len(etypes) > 0 and len(enodes) > 0:
                line_nodes = enodes[0].reshape(-1, 2)
                for ln in line_nodes:
                    pts = [coords[node_map[n]][:2] for n in ln]
                    
                    if 'fault' in name.lower():
                        fault_lines.append(pts)
                    elif 'boundary' in name.lower():
                        boundary_lines.append(pts)
    
    # Draw boundaries
    if show_boundaries and boundary_lines:
        for pts in boundary_lines:
            xs, ys = zip(*pts)
            ax.plot(xs, ys, 'g-', linewidth=1.5, zorder=5)
        legend_patches.append(Patch(facecolor='none', edgecolor='green', 
                                   linewidth=2, label='Boundaries'))
    
    # Draw faults
    if show_fault and fault_lines:
        for pts in fault_lines:
            xs, ys = zip(*pts)
            ax.plot(xs, ys, 'r-', linewidth=2.5, zorder=10)
        legend_patches.append(Patch(facecolor='none', edgecolor='red', 
                                   linewidth=2, label='Fault'))
    
    # Get mesh statistics
    _, _, _ = gmsh.model.mesh.getNodes()
    total_nodes = len(node_tags)
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(dim=2)
    total_elements = sum(len(et) for et in elem_tags)
    
    # Finalize plot
    ax.autoscale()
    ax.set_aspect('equal')
    ax.set_xlabel('X (m)', fontsize=12)
    ax.set_ylabel('Y (m)', fontsize=12)
    ax.set_title(f'{mesh_file}\n{total_nodes} nodes, {total_elements} elements', fontsize=14)
    
    # Add legend
    ax.legend(handles=legend_patches, loc='upper left', fontsize=10, 
              framealpha=0.9, edgecolor='black')
    
    # Add grid
    ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f'Saved: {output_file}')
    
    gmsh.finalize()
    
    return output_file


def main():
    parser = argparse.ArgumentParser(
        description='Plot mesh elements colored by material section',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('mesh_file', help='Input mesh file (.msh)')
    parser.add_argument('-o', '--output', help='Output image file (default: <mesh>_materials.png)')
    parser.add_argument('--no-boundaries', action='store_true', help='Hide boundary edges')
    parser.add_argument('--no-fault', action='store_true', help='Hide fault edges')
    
    args = parser.parse_args()
    
    plot_mesh(
        args.mesh_file,
        output_file=args.output,
        show_boundaries=not args.no_boundaries,
        show_fault=not args.no_fault
    )


if __name__ == '__main__':
    main()
