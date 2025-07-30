import gmsh
import numpy as np


def create_structured_nozzle_mesh(points_list, revolve_angle, cell_len, show_gui=True, 
                                  boundary_layer_thickness=0.0005, boundary_layer_ratio=1.15, boundary_layer_num=7):
    """Create structured 3D nozzle mesh with Gmsh - REVOLVE AROUND AXIS with boundary layers"""
    print("=== STEP 3: Creating Structured Mesh in Gmsh (Revolved Sector with Boundary Layers) ===")
    
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)  # Reduce verbosity (default is 5)
    gmsh.model.add("structured_nozzle_revolve")
    
    # Extract coordinates and external points
    x_coords = [point.xy_loc[0] for point in points_list]
    y_coords = [point.xy_loc[1] for point in points_list]
    external_points = [point for point in points_list if hasattr(point, 'external_points') and point.external_points]

    # Create geometry
    gmsh_points = create_points(x_coords, y_coords, external_points)
    gmsh_lines = create_lines(gmsh_points)
    gmsh_surfaces = create_surfaces(gmsh_lines)

    create_structured_mesh_with_boundary_layers(gmsh_lines, gmsh_surfaces, cell_len, 
                                               boundary_layer_thickness, boundary_layer_ratio, boundary_layer_num)
    print("Structured mesh configuration completed successfully")

    # STEP 4: Revolve 2D surfaces to create 3D volumes
    print(f"\n=== STEP 4: Revolving 2D surfaces to create 3D volumes ===")
    gmsh.model.geo.synchronize()
    
    volumes_dict = revolve_surfaces_3d(gmsh_surfaces, revolve_angle)
    print(f"Revolution completed, created {len(volumes_dict['volumes'])} volumes")
    
    # # STEP 5: Rotate domain to center the wedge
    # volumes_dict = rotate_domain(volumes_dict, revolve_angle)
    # print("Domain rotation completed")
    
    # Final synchronization before mesh generation
    gmsh.model.geo.synchronize()
    
    # Generate 3D mesh
    try:
        gmsh.model.mesh.generate(3)
        print("3D mesh generation successful")
    except Exception as e:
        print(f"3D mesh generation failed: {e}")
        print("Trying 2D then 3D mesh generation...")
        gmsh.model.mesh.generate(2)
        gmsh.model.mesh.generate(3)
        
    
    # Create boundary conditions with geometry data
    create_boundary_conditions(revolve_angle, points_list, external_points)
    
    # Get 3D mesh statistics
    nodes_3d = gmsh.model.mesh.getNodes()
    elements_3d = gmsh.model.mesh.getElements()
    
    print(f"3D Mesh Statistics:")
    print(f"  - Nodes: {len(nodes_3d[0])}")
    print(f"  - Elements: {sum(len(elem[1]) for elem in elements_3d[1:])}")
    
    # Save the mesh with OpenFOAM compatible format
    mesh_file = "structured_nozzle_3D.msh2"
     
    # Write the mesh
    gmsh.write(mesh_file)
    print(f"Mesh written to: {mesh_file}")
    print(f"Mesh format: MSH 2.2 ASCII with linear elements only (no geometric entities)")
    
    # Show GUI if requested
    if show_gui:
        print("Opening Gmsh GUI...")
        gmsh.fltk.run()
    
    gmsh.finalize()
    return "structured_nozzle_3D.msh2"
    

def create_points(x_coords, y_coords, external_points):
    '''Create points for nozzle geometry and external domain'''
    main_points = len(x_coords) - len(external_points)
    
    # Create main nozzle points (wall and centerline)
    wall_points = [gmsh.model.geo.addPoint(x_coords[i], y_coords[i], 0) for i in range(main_points)]
    centerline_points = [gmsh.model.geo.addPoint(x_coords[i], 0, 0) for i in range(main_points)]

    # Create external points
    ext_points_gmsh = [gmsh.model.geo.addPoint(ext.xy_loc[0], ext.xy_loc[1], 0) for ext in external_points]

    return {
        'wall_points': wall_points,
        'centerline': centerline_points,
        'ext_points': ext_points_gmsh
    }

def create_lines(gmsh_points):
    '''Create lines for nozzle geometry'''
    wall_pts, center_pts, ext_pts = gmsh_points['wall_points'], gmsh_points['centerline'], gmsh_points['ext_points']
    
    # Create consecutive point connections
    wall_lines = [gmsh.model.geo.addLine(wall_pts[i], wall_pts[i+1]) for i in range(len(wall_pts)-1)]
    centerline_lines = [gmsh.model.geo.addLine(center_pts[i], center_pts[i+1]) for i in range(len(center_pts)-1)]
    vertical_lines = [gmsh.model.geo.addLine(wall_pts[i], center_pts[i]) for i in range(len(wall_pts))]
    
    # External domain lines
    external_lines = [
        gmsh.model.geo.addLine(wall_pts[-1], ext_pts[0]),      # Exit to ext[0]
        gmsh.model.geo.addLine(ext_pts[0], ext_pts[1]),        # ext[0] to ext[1]
        gmsh.model.geo.addLine(ext_pts[1], center_pts[-1]),    # ext[1] to center end
        gmsh.model.geo.addLine(ext_pts[0], ext_pts[2]),        # ext[0] to ext[2]
        gmsh.model.geo.addLine(ext_pts[2], ext_pts[3]),        # ext[2] to ext[3]
        gmsh.model.geo.addLine(ext_pts[3], wall_pts[-1])       # ext[3] to exit
    ]

    return {
        'wall_lines': wall_lines,
        'centerline_lines': centerline_lines,
        'vertical_lines': vertical_lines,
        'external_lines': external_lines
    }

def create_surfaces(gmsh_lines):
    '''Create surfaces for nozzle geometry and external domain'''
    nozzle_surfaces = []

    # Create nozzle surfaces (rectangular between consecutive vertical lines)
    for i in range(len(gmsh_lines['wall_lines'])):
        try:
            curve_loop = gmsh.model.geo.addCurveLoop([
                gmsh_lines['wall_lines'][i],           # Top
                gmsh_lines['vertical_lines'][i + 1],   # Right
                -gmsh_lines['centerline_lines'][i],    # Bottom (reversed)
                -gmsh_lines['vertical_lines'][i]       # Left (reversed)
            ])
            surface = gmsh.model.geo.addPlaneSurface([curve_loop])
            nozzle_surfaces.append(surface)
        except Exception as e:
            print(f"Warning: Could not create nozzle surface {i}: {e}")

    # External surfaces
    ext_lines = gmsh_lines['external_lines']
    vert_lines = gmsh_lines['vertical_lines']
    
    ext_loop_1 = gmsh.model.geo.addCurveLoop([ext_lines[0], ext_lines[1], ext_lines[2], -vert_lines[-1]])
    ext_surface_1 = gmsh.model.geo.addPlaneSurface([ext_loop_1])

    ext_loop_2 = gmsh.model.geo.addCurveLoop([ext_lines[3], ext_lines[4], ext_lines[5], ext_lines[0]])
    ext_surface_2 = gmsh.model.geo.addPlaneSurface([ext_loop_2])
    
    return {
        'nozzle_surfaces': nozzle_surfaces,
        'external_surface_inner': ext_surface_1,
        'external_surface_outer': ext_surface_2
    }

def create_structured_mesh_with_boundary_layers(gmsh_lines, gmsh_surfaces, cell_len, 
                                               boundary_layer_thickness, boundary_layer_ratio, boundary_layer_num):
    '''Create structured mesh with inflation layers near nozzle walls'''
    
    # Synchronize the geometry first
    gmsh.model.geo.synchronize()
    
    print(f"\n=== CREATING BOUNDARY LAYER MESH ===")
    print(f"Boundary layer thickness: {boundary_layer_thickness}")
    print(f"Boundary layer growth ratio: {boundary_layer_ratio}")
    print(f"Number of boundary layer elements: {boundary_layer_num}")
    
    # Set transfinite properties for axial lines (wall and centerline)
    # FIXED: All axial lines should use individual length / cell_len for uniform x-direction cells
    for line in gmsh_lines['wall_lines'] + gmsh_lines['centerline_lines']:
        try:
            bbox = gmsh.model.getBoundingBox(1, line)
            x_min, y_min, z_min, x_max, y_max, z_max = bbox
            line_len = x_max - x_min
            
            # Each axial line should have nodes based on its own length to maintain constant cell_len
            num_nodes = max(2, int(round(line_len / cell_len)) + 1)
            print(f"Axial line {line}: length={line_len:.6f}, cell_len={cell_len:.6f}, nodes={num_nodes}")
            
            gmsh.model.geo.mesh.setTransfiniteCurve(line, num_nodes)
        except Exception as e:
            print(f"Warning: Could not set transfinite curve for line {line}: {e}")
    
    # MODIFIED: Radial lines (vertical) with boundary layer distribution
    print(f"\n=== SETTING BOUNDARY LAYER DISTRIBUTION ON VERTICAL LINES ===")
    
    # Calculate total radial length from first vertical line
    first_vertical = gmsh_lines['vertical_lines'][0]
    bbox = gmsh.model.getBoundingBox(1, first_vertical)
    x_min, y_min, z_min, x_max, y_max, z_max = bbox
    radial_length = np.sqrt((x_max - x_min)**2 + (y_max - y_min)**2 + (z_max - z_min)**2)
    
    # Calculate boundary layer node distribution
    bl_nodes, bl_progression_ratio = calculate_boundary_layer_distribution(
        radial_length, boundary_layer_thickness, boundary_layer_ratio, boundary_layer_num, cell_len)
    
    print(f"Radial length: {radial_length:.6f}")
    print(f"Total nodes in radial direction: {bl_nodes}")
    print(f"Boundary layer progression ratio: {bl_progression_ratio} (fine cells near WALL)")
    
    # Apply boundary layer distribution to ALL vertical lines
    # NOTE: Vertical lines go from WALL (y_max) to CENTERLINE (y=0)
    # We want fine cells near WALL, so we need progression > 1.0 (fine to coarse from wall to centerline)
    for i, line in enumerate(gmsh_lines['vertical_lines']):
        try:
            bbox = gmsh.model.getBoundingBox(1, line)
            x_min, y_min, z_min, x_max, y_max, z_max = bbox
            line_len = np.sqrt((x_max - x_min)**2 + (y_max - y_min)**2 + (z_max - z_min)**2)
            
            # CORRECTED: Use boundary_layer_ratio directly (not inverted)
            # This creates fine cells at the START of the line (wall) and coarse cells at END (centerline)
            gmsh.model.geo.mesh.setTransfiniteCurve(line, bl_nodes, "Progression", bl_progression_ratio)
            print(f"Vertical line {i} ({line}): length={line_len:.6f}, BL nodes={bl_nodes}, wall-to-center progression={bl_progression_ratio:.3f}")
            
        except Exception as e:
            print(f"Warning: Could not set boundary layer for vertical line {line}: {e}")
            # Fallback to uniform distribution
            try:
                uniform_nodes = max(2, int(radial_length / cell_len) + 1)
                gmsh.model.geo.mesh.setTransfiniteCurve(line, uniform_nodes)
                print(f"  -> Fallback: uniform distribution with {uniform_nodes} nodes")
            except Exception as e2:
                print(f"  -> Complete fallback failed: {e2}")

    # External lines: SAME CELL SIZE AS NOZZLE for uniform cells
    print(f"\n=== SETTING EXTERNAL LINE NODE COUNTS FOR UNIFORM CELLS ===")
    
    # Use SAME cell length as nozzle for uniform mesh throughout
    external_cell_len = cell_len  # Same as nozzle for uniform cells
    
    for i, line in enumerate(gmsh_lines['external_lines']):
        try:
            bbox = gmsh.model.getBoundingBox(1, line)
            x_min, y_min, z_min, x_max, y_max, z_max = bbox
            line_len = np.sqrt((x_max - x_min)**2 + (y_max - y_min)**2 + (z_max - z_min)**2)
            
            # Use same cell length as nozzle for uniform appearance
            num_nodes = max(2, int(round(line_len / external_cell_len)) + 1)
            print(f"External line {i} ({line}): length={line_len:.6f}, uniform_cell_len={external_cell_len:.6f}, nodes={num_nodes}")
            
            gmsh.model.geo.mesh.setTransfiniteCurve(line, num_nodes)
        except Exception as e:
            print(f"Warning: Could not set transfinite curve for external line {line}: {e}")
    
    # Set same mesh size for external points
    try:
        # Get external points from geometry
        all_points = gmsh.model.getEntities(0)  # Get all points
        
        # Set mesh size for points that are likely external (high x or y coordinates)
        for dim, point_tag in all_points:
            bbox = gmsh.model.getBoundingBox(dim, point_tag)
            x_min, y_min, z_min, x_max, y_max, z_max = bbox
            
            # If point is far from origin, it's likely external - set same mesh size as nozzle
            if abs(x_min) > 0.1 or abs(y_min) > 0.1:
                gmsh.model.geo.mesh.setSize([(dim, point_tag)], external_cell_len)
    except Exception as e:
        print(f"Warning: Could not set mesh sizes for external points: {e}")
    
    # Set structured mesh for surfaces with boundary layer consideration
    print(f"\n=== SETTING TRANSFINITE SURFACES FOR STRUCTURED MESH WITH BOUNDARY LAYERS ===")
    
    # Process nozzle surfaces (these should work correctly with boundary layers)
    for i, surface in enumerate(gmsh_surfaces['nozzle_surfaces']):
        if surface:
            try:
                # Set transfinite surface for true structured mesh
                gmsh.model.geo.mesh.setTransfiniteSurface(surface)
                # Set recombine to get quad elements
                gmsh.model.geo.mesh.setRecombine(2, surface)
                print(f"Nozzle surface {i} ({surface}): Set as transfinite structured surface with boundary layers")
            except Exception as e:
                print(f"Warning: Could not set transfinite nozzle surface {surface}: {e}")
    
    # Handle external surfaces - ONLY RECOMBINE (uniform cell size but no transfinite)
    external_surfaces = [gmsh_surfaces['external_surface_inner'], gmsh_surfaces['external_surface_outer']]
    for i, surface in enumerate(external_surfaces):
        if surface:
            try:
                # Use structured meshing algorithm + recombine for uniform appearance
                gmsh.model.geo.mesh.setAlgorithm(2, surface, 8)  # Frontal-Delaunay for quads
                gmsh.model.geo.mesh.setRecombine(2, surface)
                print(f"External surface {i} ({surface}): Set with structured algorithm + recombine (uniform cells)")
            except Exception as e:
                print(f"Warning: External surface {i} algorithm failed: {e}")
                try:
                    # Fallback: just recombine
                    gmsh.model.geo.mesh.setRecombine(2, surface)
                    print(f"External surface {i}: Basic recombine only")
                except Exception as e2:
                    print(f"External surface {i} will use default meshing: {e2}")
    
    return gmsh_surfaces

def calculate_boundary_layer_distribution(radial_length, boundary_layer_thickness, boundary_layer_ratio, boundary_layer_num, cell_len):
    '''Calculate boundary layer node distribution parameters'''
    
    # Validate inputs
    if boundary_layer_thickness >= radial_length:
        print(f"Warning: Boundary layer thickness ({boundary_layer_thickness}) >= radial length ({radial_length})")
        print("  -> Reducing boundary layer thickness to 60% of radial length")
        boundary_layer_thickness = 0.6 * radial_length
    
    if boundary_layer_ratio <= 1.0:
        print(f"Warning: Boundary layer ratio must be > 1.0, got {boundary_layer_ratio}")
        boundary_layer_ratio = 1.15  # Less aggressive than 1.2
        print(f"  -> Using default ratio: {boundary_layer_ratio}")
    
    # Make boundary layer less dense - reduce number of BL cells if too many
    if boundary_layer_num > radial_length / (cell_len * 0.1):  # Don't use more than 10% of available cells for BL
        boundary_layer_num = max(5, int(radial_length / (cell_len * 0.2)))  # Use 20% of cells for BL
        print(f"  -> Reducing BL cells to {boundary_layer_num} for less dense mesh")
    
    # Calculate boundary layer region parameters
    # First layer thickness (closest to wall) - use a more conservative approach
    first_layer_thickness = boundary_layer_thickness / sum([boundary_layer_ratio**i for i in range(boundary_layer_num)])
    
    # Total thickness of boundary layer region
    bl_total_thickness = first_layer_thickness * (boundary_layer_ratio**boundary_layer_num - 1) / (boundary_layer_ratio - 1)
    
    # Check if boundary layer fits in radial direction
    if bl_total_thickness > 0.7 * radial_length:  # Use max 70% of radial space for BL
        # Adjust boundary layer to fit
        scale_factor = 0.6 * radial_length / bl_total_thickness
        boundary_layer_thickness *= scale_factor
        bl_total_thickness *= scale_factor
        first_layer_thickness *= scale_factor
        print(f"Warning: Boundary layer rescaled by factor {scale_factor:.3f} to fit geometry")
    
    # Remaining region after boundary layer
    outer_region_length = radial_length - bl_total_thickness
    
    # Nodes in outer region (uniform distribution) - ensure we have reasonable outer region
    outer_nodes = max(3, int(outer_region_length / cell_len))  # At least 3 nodes in outer region
    
    # Total nodes = boundary layer nodes + outer region nodes
    total_nodes = boundary_layer_num + outer_nodes + 1  # +1 for the wall node
    
    # Calculate y+ estimation for reference
    # Assuming typical nozzle flow conditions
    velocity_estimate = 100  # m/s (rough estimate)
    density_estimate = 1.0   # kg/m³
    viscosity_estimate = 1.8e-5  # kg/(m·s) for air
    wall_shear_stress = 0.5 * density_estimate * velocity_estimate**2 * 0.005  # Rough Cf estimate
    u_tau = np.sqrt(wall_shear_stress / density_estimate)
    y_plus_first_cell = first_layer_thickness * u_tau * density_estimate / viscosity_estimate
    
    print(f"  First layer thickness: {first_layer_thickness:.6f} m")
    print(f"  BL total thickness: {bl_total_thickness:.6f} m ({100*bl_total_thickness/radial_length:.1f}% of radial)")
    print(f"  Outer region length: {outer_region_length:.6f} m ({100*outer_region_length/radial_length:.1f}% of radial)")
    print(f"  BL cells: {boundary_layer_num}, Outer cells: {outer_nodes}")
    print(f"  Estimated y+ of first cell: {y_plus_first_cell:.1f} (target: <1 for LES, <30 for RANS)")
    
    return total_nodes, boundary_layer_ratio  # Return the actual ratio, not inverted

def revolve_surfaces_3d(gmsh_surfaces, revolve_angle=1.0):
    """Revolve 2D surfaces around x-axis to create 3D volumes"""
    
    angle_rad = revolve_angle * np.pi / 180.0
    all_surfaces = (gmsh_surfaces['nozzle_surfaces'] + 
                   [gmsh_surfaces['external_surface_inner'], gmsh_surfaces['external_surface_outer']])
    
    volumes = []
    for surface in all_surfaces:
        if surface:
            try:
                revolved_entities = gmsh.model.geo.revolve(
                    [(2, surface)], 0, 0, 0, 1, 0, 0, angle_rad, [1], [1.0], True)
                
                for dim, tag in revolved_entities:
                    if dim == 3:
                        volumes.append(tag)
                        gmsh.model.geo.mesh.setTransfiniteVolume(tag)
            except Exception as e:
                print(f"Warning: Could not revolve surface {surface}: {e}")
    
    return {'volumes': volumes}

def rotate_domain(volumes_dict, revolve_angle=1.0):
    """Rotate the entire 3D domain by -revolve_angle/2 around x-axis to center the wedge"""
    
    volumes = volumes_dict.get('volumes', [])
    half_angle_rad = -(revolve_angle / 2.0) * np.pi / 180.0
    
    # Get all entities, fallback to volumes only
    try:
        all_entities = sum([gmsh.model.geo.getEntities(dim) for dim in [0,1,2,3]], [])
        print(f"Rotating {len(all_entities)} entities by -{revolve_angle/2}°")
    except Exception as e:
        all_entities = [(3, vol) for vol in volumes]
        print(f"Warning: Volume fallback rotation: {e}")
    
    # Apply rotation
    try:
        if all_entities:
            gmsh.model.geo.rotate(all_entities, 0, 0, 0, 1, 0, 0, half_angle_rad)
    except Exception as e:
        print(f"Error: Rotation failed: {e}")
    
    return volumes_dict

def create_boundary_conditions(revolve_angle=1.0, points_list=None, external_points=None):

    print("=== Creating boundary conditions ===")

    # Get all entities in the model
    all_surfaces = gmsh.model.getEntities(2)  # 2D entities (surfaces)
    all_volumes = gmsh.model.getEntities(3)   # 3D entities (volumes)
    all_lines = gmsh.model.getEntities(1)     # 1D entities (lines) for wedge/axis

    # Analyze ALL surfaces using geometric properties
    surface_data_dict = get_all_surface_data(all_surfaces)
    
    # Classify surfaces into boundary conditions
    frontandback_surfaces = []
    inlet_surfaces = []
    outlet_surfaces = []
    wall_surfaces = []
    
    # Classify each surface
    for surface_tag, surface_data in surface_data_dict.items():
        print(f"\nAnalyzing surface {surface_tag} for boundary conditions:")
        if analyze_frontandback_boundary(surface_data, revolve_angle):
            frontandback_surfaces.append(surface_tag)
        elif analyze_inlet_boundary(surface_data, points_list):
            inlet_surfaces.append(surface_tag)
        elif analyze_outlet_boundary(surface_data, points_list) == True:
            outlet_surfaces.append(surface_tag)
        elif analyze_outlet_boundary(surface_data, points_list) == "skip":
            print(f"    -> Skipping surface {surface_tag} (domain point, not boundary)")
            continue  # Skip this surface entirely
        elif analyze_internal_boundary(surface_data):
            continue
        else:
             wall_surfaces.append(surface_tag)

    # Get wedge lines (centerline)
    wedge_lines = analyze_all_lines_for_wedge(all_lines)
    
    boundary_conditions = []

    # Store surfaces for each boundary condition
    boundary_conditions.append({
        'name': 'frontAndBack',
        'key': 'front_and_back',
        'dimension': 2,
        'entity_type': 'surfaces',
        'entities': frontandback_surfaces
    })
    boundary_conditions.append({
        'name': 'inlet',
        'key': 'inlet',
        'dimension': 2,
        'entity_type': 'surfaces',
        'entities': inlet_surfaces
    })
    boundary_conditions.append({
        'name': 'outlet',
        'key': 'outlet',
        'dimension': 2,
        'entity_type': 'surfaces',
        'entities': outlet_surfaces
    })
    boundary_conditions.append({
        'name': 'wall',
        'key': 'wall',
        'dimension': 2,
        'entity_type': 'surfaces',
        'entities': wall_surfaces
    })
    boundary_conditions.append({
        'name': 'wedge',
        'key': 'wedge',
        'dimension': 1,
        'entity_type': 'lines',
        'entities': wedge_lines
    })
    boundary_conditions.append({
        'name': 'fluid',
        'key': 'fluid',
        'dimension': 3,
        'entity_type': 'volumes',
        'entities': [tag for dim, tag in all_volumes]
    })

    
    # Create physical groups for OpenFOAM
    physical_groups = create_physical_groups_for_openfoam(boundary_conditions)
        
    return physical_groups


    
def get_all_surface_data(all_surfaces):
    """
    Analyze ALL surfaces using geometric properties instead of mesh nodes
    This works for both meshed and unmeshed surfaces
    Returns a dictionary with surface tag as key and surface data as value
    """
    print("\n=== GEOMETRIC SURFACE ANALYSIS (ALL SURFACES) ===")

    surface_data_dict = {}

    # Analyze each surface
    for dim, tag in all_surfaces:
        try:
            # Get geometric bounding box of the surface
            bbox = gmsh.model.getBoundingBox(dim, tag)
            x_min, y_min, z_min, x_max, y_max, z_max = bbox
            
            size_xyz = { 'x':x_max - x_min, 'y':y_max - y_min, 'z':z_max - z_min}
            print(f"\nSurface {tag} geometric analysis:")
            print(f"  Size: ({size_xyz['x']:.6f}, {size_xyz['y']:.6f}, {size_xyz['z']:.6f})")
            
            # Store surface data
            surface_data_dict[tag] = {
                'tag': tag,
                'size': (size_xyz['x'], size_xyz['y'], size_xyz['z']),
                'bbox': bbox
            }
                
        except Exception as e:
            print(f"  Warning: Could not analyze surface {tag}: {e}")
            continue
        
    return surface_data_dict    


def analyze_frontandback_boundary(surface_data, revolve_angle):
    """Analyze if surface is a frontAndBack boundary based on normal vector angle
    For revolved geometry, front/back surfaces are the original 2D surfaces at z = ±half_angle"""
    size_x, size_y, size_z = surface_data['size']
    
    if abs(np.degrees(np.arctan(size_z / size_y)) - revolve_angle ) < 0.0001:  # Avoid division by zero
        print(f"    -> FRONT/BACK surface detected ")
        return True
    
    if abs(np.degrees(np.arctan(size_z / size_y))) < 0.0001:  # Avoid division by zero
        print(f"    -> FRONT/BACK surface detected ")
        return True
    
    return False

def analyze_inlet_boundary(surface_data, points_list):
    """Analyze if surface is an inlet boundary (at minimum x position)"""   
    x_min = np.min([point.xy_loc[0] for point in points_list])
    tolerance = 1e-6
    if abs(surface_data['bbox'][0] - x_min) < tolerance and surface_data['size'][0] < tolerance:
        print(f"    -> INLET surface detected (x_min={surface_data['bbox'][0]:.6f}, domain_x_min={x_min:.6f})")
        return True
    return False

def analyze_internal_boundary(surface_data):
    """
    Analyze if surface is a wall boundary using nozzle wall points
    """
    size_x, size_y, size_z = surface_data['size']
    print('wall =',abs(np.degrees(np.arctan(size_z / size_y))))
    if abs(np.degrees(np.arctan(size_z / size_y)))   < 1.01 and abs(np.degrees(np.arctan(size_z / size_y))) >0.998:  # Avoid division by zero
        print(f"    -> internal surface ")
        return True
    
    return False

def analyze_outlet_boundary(surface_data, points_list):
    """Analyze if surface is an outlet boundary"""
    x_max_domain = np.max([point.xy_loc[0] for point in points_list])
    
    # Find the nozzle wall point with maximum x-coordinate
    nozzle_wall_points = [point for point in points_list if hasattr(point, 'nozzle_walls') and point.nozzle_walls]
    nozzle_end_point = max(nozzle_wall_points, key=lambda p: p.xy_loc[0])
    y_nozzle_end = nozzle_end_point.xy_loc[1]
    x_nozzle_end = nozzle_end_point.xy_loc[0]

    
    # Get y-coordinates of external points only
    external_y_coords = [point.xy_loc[1] for point in points_list 
                        if hasattr(point, 'external_points') and point.external_points]
    y_max_domain = np.max(external_y_coords) if external_y_coords else 0.0

    x_min, y_min, z_min, x_max, y_max, z_max = surface_data['bbox']
    
    tolerance = 1e-3

    if x_max < x_max_domain + tolerance and x_max > x_max_domain - tolerance:
        print(f"    -> OUTLET surface detected ")
        if (y_nozzle_end < y_min + tolerance and y_nozzle_end > y_min - tolerance) and \
           (x_nozzle_end < x_min + tolerance and x_nozzle_end > x_min - tolerance):
            return "skip"  # Special return value to skip this surface entirely
        else:
            return True
    
    if y_max < y_max_domain + tolerance and y_max > y_max_domain - tolerance:
        print(f"    -> OUTLET surface detected ")
        return True
    
    
    return False

def analyze_all_lines_for_wedge(all_lines):
    """
    Analyze ALL lines to find centerline lines for wedge boundary condition
    Wedge BC should be applied to lines (1D entities) that lie on centerline (y ≈ 0)
    """
    print("\n=== LINE ANALYSIS FOR WEDGE BOUNDARY CONDITION ===")
    
    wedge_lines = []  # Lines on centerline (y ≈ 0) for wedge BC
    
    for dim, tag in all_lines:
        # Get geometric bounding box of the line
        bbox = gmsh.model.getBoundingBox(dim, tag)
        x_min, y_min, z_min, x_max, y_max, z_max = bbox
                    
        size_y = y_max - y_min
        
        # Check if line is on centerline (y ≈ 0)
        tolerance = 1e-6
        
        if abs(y_min) < tolerance and abs(y_max) < tolerance and size_y < tolerance:
            wedge_lines.append(tag)
            # print(f"  -> WEDGE LINE (centerline)")

    return wedge_lines

def create_physical_groups_for_openfoam(boundary_conditions):
    """
    Create physical groups for OpenFOAM boundary conditions
    """
    print("\n=== CREATING PHYSICAL GROUPS FOR OPENFOAM ===")
    
    # Create physical groups for boundaries
    physical_groups = {}
    
    # Create physical groups for all boundary conditions
    for bc in boundary_conditions:
        if bc['entities']:  # Check if there are entities for this boundary condition
            entities = bc['entities']
            group = gmsh.model.addPhysicalGroup(bc['dimension'], entities)
            gmsh.model.setPhysicalName(bc['dimension'], group, bc['name'])
            physical_groups[bc['name']] = group
            print(f"Created '{bc['name']}' physical group with {len(entities)} {bc['entity_type']}")
        else:
            print(f"Warning: No {bc['entity_type']} found for '{bc['name']}' boundary condition!")
            if bc['name'] == 'wedge':
                print("This might cause issues with OpenFOAM wedge boundary condition.")
                print("Wedge BC requires centerline lines (y ≈ 0) for axisymmetric simulations.")
    
    print(f"Total physical groups created: {len(physical_groups)}")
    return physical_groups

