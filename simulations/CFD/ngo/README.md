# CFD Nozzle Analysis

A comprehensive CFD (Computational Fluid Dynamics) toolkit for nozzle design and analysis using Gmsh and OpenFOAM.

## Features

- **Structured Mesh Generation**: Create high-quality structured 3D meshes for nozzle geometries using Gmsh
- **Boundary Condition Management**: Automatic detection and assignment of boundary conditions for OpenFOAM
- **OpenFOAM Integration**: Complete workflow for CFD simulations with OpenFOAM
- **Nozzle Geometry Processing**: Tools for processing and analyzing nozzle geometries

## Main Components

### Core Scripts

- `nozzleGMSH.py` - Main mesh generation script with advanced boundary condition detection
- `cfdmain.py` - Main CFD workflow coordinator
- `openfoam_preprocessing.py` - OpenFOAM case preparation and preprocessing
- `complete_workflow.py` - End-to-end automation script

### Key Features in nozzleGMSH.py

- **Structured Mesh Creation**: Generates high-quality hexahedral meshes
- **Revolve Around Axis**: Creates 3D meshes by revolving 2D profiles
- **Boundary Detection**: Automatic classification of surfaces into:
  - Inlet/Outlet boundaries
  - Wall boundaries
  - Front/Back surfaces (for wedge simulations)
  - Centerline (wedge boundary condition)
- **OpenFOAM Compatibility**: Proper physical group creation for OpenFOAM

## Usage

### Basic Mesh Generation

```python
from nozzleGMSH import create_structured_nozzle_mesh

# Generate mesh from nozzle points with boundary layers
mesh_file = create_structured_nozzle_mesh(
    points_list=nozzle_points,
    revolve_angle=1.0,  # 1 degree wedge
    cell_len=0.001,     # Base cell length
    show_gui=True,
    boundary_layer_thickness=0.0005,  # 0.5mm boundary layer (less dense)
    boundary_layer_ratio=1.15,        # Gentler growth ratio
    boundary_layer_num=7               # Fewer BL cells
)
```

### Boundary Conditions

The mesh generator automatically creates these boundary conditions:
- `inlet` - Chamber inlet
- `outlet` - Nozzle outlet  
- `wall` - Nozzle walls
- `frontAndBack` - Wedge front/back surfaces
- `wedge` - Centerline for axisymmetric simulations
- `fluid` - Fluid volume

## Requirements

- Python 3.7+
- Gmsh Python API
- NumPy
- OpenFOAM (for CFD simulations)

## Test Cases

The `test/` directory contains complete OpenFOAM test cases with:
- Boundary conditions setup
- Solver configurations
- Post-processing scripts

## License

This project is part of a rocket engine design and simulation suite.
