#!/usr/bin/env python3
"""
CFD Main Script - Complete Nozzle Geometry and Mesh Generation Workflow
This script:
1. Generates nozzle geometry points
2. Plots the nozzle contour
3. Creates a structured 3D mesh in Gmsh
4. Displays the mesh in Gmsh GUI
5. Provides OpenFOAM conversion instructions
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
import sys
import os
from . import nozzleGMSH
from . import preprocessingNGO 

@dataclass
class InitialisePoint:
    idx: int
    xy_loc: list[float] = field(default_factory=lambda: [0.0, 0.0])
    external_points: bool = False  # Flag external points
    nozzle_walls: bool = False  # Flag nozzle walls
    chamber_points: bool = False  # Flag chamber points

def generate_nozzle_geometry(y_combustion_chamber, y_exit, y_throat, num_of_points_nozzle):
    """Generate nozzle geometry points including combustion chamber and downstream domain"""
    print("=== STEP 1: Generating Nozzle Geometry with Combustion Chamber and Downstream Domain ===")
    print(f"Parameters: y_combustion_chamber={y_combustion_chamber}, y_exit={y_exit}, y_throat={y_throat}, num_of_points_nozzle={num_of_points_nozzle}")
    # Nozzle parameters
    num_of_points = int(num_of_points_nozzle * 0.9) # Number of points in the nozzle
    angle_inlet = 45    # Inlet angle (degrees)
    angle_outlet = 16   # Outlet angle (degrees)
    
    # Combustion chamber parameters
    num_chamber_points = int(num_of_points_nozzle * 0.1)  # Fewer points for shorter chamber
    
    # Geometry setup
    ymid = y_throat  # Throat half-height
    ystart = y_combustion_chamber  # Starting height (chamber height)
    yend = y_exit    # Outlet half-height
    
    # Combustion chamber starts at negative x
    xchamber_start = -y_combustion_chamber
    xchamber_end = 0.0  # Chamber ends at nozzle inlet
    
    # Nozzle geometry
    xnozzle_start = xchamber_end
    xmid = xnozzle_start + (ystart - ymid) / np.tan(np.radians(angle_inlet))
    xnozzle_end = xmid + (yend - ymid) / np.tan(np.radians(angle_outlet))
    
    # Downstream domain
    xdownstream_end = (xnozzle_end + (xnozzle_end- xnozzle_start)) * 2
    y_offset = yend * 1  # Offset by twice the nozzle exit height

    points_list = []

    # 1. Generate combustion chamber points
    deltax_chamber = abs(xchamber_start / num_chamber_points)
    for i in range(num_chamber_points):
        point = InitialisePoint(idx=len(points_list))
        x = xchamber_start + i * deltax_chamber  # Fixed: use i instead of num_chamber_points
        y = ystart  # Constant height in chamber
        point.xy_loc = [x, y]
        point.nozzle_walls = True
        points_list.append(point)
        point.chamber_points = True

    # 2. Generate nozzle points
    deltax_nozzle = (xnozzle_end - xnozzle_start) / num_of_points
    for i in range(num_of_points + 1):  # +1 to include end point
        point = InitialisePoint(idx=len(points_list))
        
        if i == 0:  # First point (nozzle inlet)
            x, y = xnozzle_start, ystart
        elif i == num_of_points:  # Last point (nozzle outlet)
            x, y = xnozzle_end, yend
        else:
            x = xnozzle_start + i * deltax_nozzle
            if x <= xmid:  # Converging section
                y = ystart - (x - xnozzle_start) * np.tan(np.radians(angle_inlet))
                point.nozzle_converging_section = True
            else:  # Diverging section
                y = ymid + (x - xmid) * np.tan(np.radians(angle_outlet))
                point.nozzle_diverging_section = True
        
        point.xy_loc = [x, y]
        point.nozzle_walls = True  # Flag all nozzle points as nozzle walls
        points_list.append(point)

    # 1. external point
    point = InitialisePoint(idx=len(points_list))
    point.xy_loc = [xdownstream_end , yend]
    point.external_points = True  # Flag as external point
    points_list.append(point)

    # 2. external point
    point_corner = InitialisePoint(idx=len(points_list))
    point_corner.xy_loc = [xdownstream_end, 0.0]  # Same X as downstream end, but on centerline
    point_corner.external_points = True  # Flag as external point
    points_list.append(point_corner)

    # 3. external point
    point_ext2 = InitialisePoint(idx=len(points_list))
    point_ext2.xy_loc = [xdownstream_end, yend + y_offset]
    point_ext2.external_points = True  # Flag as external point
    points_list.append(point_ext2)

    # 4. external point
    point_ext1 = InitialisePoint(idx=len(points_list))
    point_ext1.xy_loc = [xnozzle_end, yend + y_offset]
    point_ext1.external_points = True  # Flag as external point
    points_list.append(point_ext1)

    nozzle_wall_points = [p for p in points_list if p.nozzle_walls]
    throat_point = min(nozzle_wall_points, key=lambda p: p.xy_loc[1])
    x_throat = throat_point.xy_loc[0]

    # 5. external point
    point_ext5 = InitialisePoint(idx=len(points_list))
    point_ext5.xy_loc = [x_throat, yend + y_offset]
    point_ext5.external_points = True  # Flag as external point
    points_list.append(point_ext5)

    # 6. external point
    point_ext6 = InitialisePoint(idx=len(points_list))
    point_ext6.xy_loc = [x_throat, yend]
    point_ext6.external_points = True  # Flag as external point
    points_list.append(point_ext6)

    print(x_throat, yend)
    print(x_throat, yend + y_offset)

    print(f"Generated {len(points_list)} geometry points")
    print(f"  - Chamber: {num_chamber_points} points")
    print(f"  - Nozzle: {num_of_points + 1} points") 
    print(f"  - External boundary: 3 points (Y-offset: {y_offset:.1f})")
    
    # Count flagged points
    external_count = sum(1 for p in points_list if p.external_points)
    nozzle_wall_count = sum(1 for p in points_list if p.nozzle_walls)
    print(f"  - External points flagged: {external_count}")
    print(f"  - Nozzle wall points flagged: {nozzle_wall_count}\n")

    return points_list


def calculate_cell_length(points_list, num_intermediate_points):
    """
    Calculate cell lengths for combustion chamber and nozzle wall sections
    
    Args:
        points_list: List of InitialisePoint objects containing geometry points
    
    Returns:
        tuple: (len_chamber_points, len_wall_nozzle_points) - average cell lengths
    """
    
    # Separate chamber points and nozzle wall points
    chamber_points = [point for point in points_list if hasattr(point, 'chamber_points') and point.chamber_points]
    nozzle_wall_points = [point for point in points_list if point.nozzle_walls and not (hasattr(point, 'chamber_points') and point.chamber_points)]
    
    # Calculate chamber cell length
    len_chamber_points = 0.0
    if len(chamber_points) > 1:
        chamber_distances = []
        for i in range(len(chamber_points) - 1):
            x1, y1 = chamber_points[i].xy_loc
            x2, y2 = chamber_points[i + 1].xy_loc
            distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
            chamber_distances.append(distance)
        len_chamber_points = np.mean(chamber_distances) if chamber_distances else 0.0
    
    # Calculate nozzle wall cell length
    len_wall_nozzle_points = 0.0
    if len(nozzle_wall_points) > 1:
        nozzle_distances = []
        for i in range(len(nozzle_wall_points) - 1):
            x1, y1 = nozzle_wall_points[i].xy_loc
            x2, y2 = nozzle_wall_points[i + 1].xy_loc
            distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
            nozzle_distances.append(distance)
        len_wall_nozzle_points = np.mean(nozzle_distances) if nozzle_distances else 0.0
    
    print(f"Cell length calculations:")
    print(f"  - Chamber points: {len(chamber_points)}, Average cell length: {len_chamber_points:.6f}")
    print(f"  - Nozzle wall points: {len(nozzle_wall_points)}, Average cell length: {len_wall_nozzle_points:.6f}")
    
    return min(len_chamber_points, len_wall_nozzle_points)/ num_intermediate_points

def start_optimisation(y_combustion_chamber, y_exit, y_throat, num_of_points_nozzle):
    """Main execution function"""
    print("CFD NOZZLE GEOMETRY AND MESH GENERATION WORKFLOW")
    print("=" * 55)
    num_intermediate_points = 3
    try:
        # Step 1: Generate geometry
        points_list = generate_nozzle_geometry(y_combustion_chamber, y_exit, y_throat, num_of_points_nozzle)
        
        cell_len = calculate_cell_length(points_list, num_intermediate_points)
        # print(len_chamber_points, len_wall_nozzle_points)

        # Step 3 & 4: Create mesh with boundary layers for accurate wall flow simulation
        # Boundary layer parameters for high-Reynolds number flow (less dense)
        boundary_layer_thickness = 0.0001  # 0.3mm boundary layer thickness (reduced further)
        boundary_layer_ratio = 1.015         # Very gentle growth ratio (10% per cell)
        boundary_layer_num = 5             # Even fewer cells in boundary layer
        
        mesh_file = nozzleGMSH.create_structured_nozzle_mesh(
            points_list, revolve_angle=1.0, cell_len=cell_len, show_gui=False,
            boundary_layer_thickness=boundary_layer_thickness,
            boundary_layer_ratio=boundary_layer_ratio, 
            boundary_layer_num=boundary_layer_num)
        
        preprocessingNGO.move_and_generate_mesh()

    except Exception as e:
        print(f"\nERROR: {e}")
        sys.exit(1)


