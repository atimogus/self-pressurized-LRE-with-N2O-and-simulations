#!/usr/bin/env python3
"""
OpenFOAM Preprocessing Script for Nozzle Flow Simulation
This script:
1. Converts Gmsh mesh to OpenFOAM format using gmshToFoam
2. Sets up a basic sonicFoam case structure
3. Creates necessary dictionaries for laminar flow simulation
4. Includes changeDictionaryDict for boundary condition modifications
"""

import os
import sys
import subprocess
import shutil
from pathlib import Path
import foamlib

class OpenFOAMNozzlePreprocessor:
    def __init__(self, mesh_file="structured_nozzle_3D.msh", case_name="nozzle_case"):
        """
        Initialize the OpenFOAM preprocessor
        
        Args:
            mesh_file: Path to the Gmsh mesh file
            case_name: Name of the OpenFOAM case directory
        """
        self.mesh_file = mesh_file
        self.case_name = case_name
        self.case_dir = Path(case_name)
        self.openfoam_version = "2412"
        
        print(f"OpenFOAM Nozzle Preprocessor initialized")
        print(f"  - Mesh file: {self.mesh_file}")
        print(f"  - Case name: {self.case_name}")
        print(f"  - OpenFOAM version: {self.openfoam_version}")
    
    def setup_case_directory(self):
        """Create the basic OpenFOAM case directory structure"""
        print("\n=== Setting up OpenFOAM case directory structure ===")
        
        # Remove existing case directory if it exists
        if self.case_dir.exists():
            shutil.rmtree(self.case_dir)
            print(f"Removed existing case directory: {self.case_dir}")
        
        # Create main directories
        directories = [
            "0",
            "constant",
            "system"
        ]
        
        for directory in directories:
            (self.case_dir / directory).mkdir(parents=True, exist_ok=True)
            print(f"Created directory: {self.case_dir}/{directory}")
    
    def convert_mesh_to_openfoam(self):
        """Convert Gmsh mesh to OpenFOAM format using gmshToFoam"""
        print("\n=== Converting Gmsh mesh to OpenFOAM format ===")
        
        if not Path(self.mesh_file).exists():
            raise FileNotFoundError(f"Mesh file not found: {self.mesh_file}")
        
        # Copy mesh file to case directory
        mesh_dest = self.case_dir / self.mesh_file
        shutil.copy2(self.mesh_file, mesh_dest)
        print(f"Copied mesh file to: {mesh_dest}")
        
        # Run gmshToFoam
        cmd = [
            "gmshToFoam",
            "-case", str(self.case_dir),
            str(mesh_dest)
        ]
        
        print(f"Running: {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            print("gmshToFoam completed successfully")
            print("Output:", result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"Error running gmshToFoam: {e}")
            print("Error output:", e.stderr)
            raise
    
    def create_control_dict(self):
        """Create the controlDict file for sonicFoam simulation"""
        print("\n=== Creating controlDict ===")
        
        control_dict_content = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     sonicFoam;

startFrom       latestTime;

startTime       0;

stopAt          endTime;

endTime         0.1;

deltaT          1e-6;

writeControl    adjustableRunTime;

writeInterval   0.01;

purgeWrite      0;

writeFormat     ascii;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

adjustTimeStep  yes;

maxCo           0.5;

maxDeltaT       1e-4;

// ************************************************************************* //
"""
        
        control_dict_path = self.case_dir / "system" / "controlDict"
        with open(control_dict_path, 'w') as f:
            f.write(control_dict_content)
        print(f"Created: {control_dict_path}")
    
    def create_fv_schemes(self):
        """Create the fvSchemes file"""
        print("\n=== Creating fvSchemes ===")
        
        fv_schemes_content = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ddtSchemes
{
    default         Euler;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
    div(phi,U)      Gauss upwind;
    div(phid,p)     Gauss upwind;
    div(phi,K)      Gauss upwind;
    div(phi,h)      Gauss upwind;
    div(phi,k)      Gauss upwind;
    div(phi,epsilon) Gauss upwind;
    div(phi,R)      Gauss upwind;
    div(R)          Gauss linear;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

// ************************************************************************* //
"""
        
        fv_schemes_path = self.case_dir / "system" / "fvSchemes"
        with open(fv_schemes_path, 'w') as f:
            f.write(fv_schemes_content)
        print(f"Created: {fv_schemes_path}")
    
    def create_fv_solution(self):
        """Create the fvSolution file"""
        print("\n=== Creating fvSolution ===")
        
        fv_solution_content = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solvers
{
    p
    {
        solver          GAMG;
        tolerance       1e-6;
        relTol          0.01;
        smoother        GaussSeidel;
    }

    pFinal
    {
        $p;
        relTol          0;
    }

    "(rho|U|h|k|epsilon)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-6;
        relTol          0.1;
    }

    "(rho|U|h|k|epsilon)Final"
    {
        $U;
        relTol          0;
    }
}

PIMPLE
{
    momentumPredictor   yes;
    nOuterCorrectors    50;
    nCorrectors         1;
    nNonOrthogonalCorrectors 0;

    residualControl
    {
        p               1e-5;
        U               1e-5;
        h               1e-5;
        
        // possibly check turbulence fields
        "(k|epsilon|omega)" 1e-5;
    }
}

relaxationFactors
{
    equations
    {
        U               0.9;
        h               0.9;
        p               0.9;
    }
}

// ************************************************************************* //
"""
        
        fv_solution_path = self.case_dir / "system" / "fvSolution"
        with open(fv_solution_path, 'w') as f:
            f.write(fv_solution_content)
        print(f"Created: {fv_solution_path}")
    
    def create_change_dictionary_dict(self):
        """Create the changeDictionaryDict file for boundary modifications"""
        print("\n=== Creating changeDictionaryDict ===")
        
        change_dict_content = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      changeDictionaryDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

boundary
{
    // Change specific boundaries from wall to patch types
    inlet
    {
        type            patch;
    }
    
    outlet
    {
        type            patch;
    }
    
    // Example: Change walls to appropriate types
    // Uncomment and modify as needed for your specific mesh
    /*
    nozzle_wall
    {
        type            wall;
    }
    
    chamber_wall
    {
        type            wall;
    }
    
    symmetry_plane
    {
        type            symmetryPlane;
    }
    
    wedge_front
    {
        type            wedge;
    }
    
    wedge_back
    {
        type            wedge;
    }
    */
}

// You can also modify field values here if needed
/*
U
{
    internalField   uniform (0 0 0);

    boundaryField
    {
        inlet
        {
            type            fixedValue;
            value           uniform (100 0 0);  // Modify as needed
        }
        
        outlet
        {
            type            zeroGradient;
        }
        
        ".*wall.*"
        {
            type            noSlip;
        }
    }
}

T
{
    internalField   uniform 300;

    boundaryField
    {
        inlet
        {
            type            fixedValue;
            value           uniform 600;  // Modify as needed
        }
        
        outlet
        {
            type            zeroGradient;
        }
        
        ".*wall.*"
        {
            type            zeroGradient;
        }
    }
}

p
{
    internalField   uniform 101325;

    boundaryField
    {
        inlet
        {
            type            zeroGradient;
        }
        
        outlet
        {
            type            fixedValue;
            value           uniform 101325;  // Modify as needed
        }
        
        ".*wall.*"
        {
            type            zeroGradient;
        }
    }
}
*/

// ************************************************************************* //
"""
        
        change_dict_path = self.case_dir / "system" / "changeDictionaryDict"
        with open(change_dict_path, 'w') as f:
            f.write(change_dict_content)
        print(f"Created: {change_dict_path}")
    
    def create_thermophysical_properties(self):
        """Create basic thermophysicalProperties file"""
        print("\n=== Creating thermophysicalProperties ===")
        
        thermo_props_content = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      thermophysicalProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

thermoType
{
    type            hePsiThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleEnthalpy;
}

mixture
{
    specie
    {
        molWeight       28.9;  // Air molecular weight
    }
    thermodynamics
    {
        Cp              1005;  // Specific heat at constant pressure [J/kg/K]
        Hf              0;     // Heat of formation [J/kg]
    }
    transport
    {
        mu              1.84e-05;  // Dynamic viscosity [Pa.s]
        Pr              0.7;       // Prandtl number
    }
}

// ************************************************************************* //
"""
        
        thermo_props_path = self.case_dir / "constant" / "thermophysicalProperties"
        with open(thermo_props_path, 'w') as f:
            f.write(thermo_props_content)
        print(f"Created: {thermo_props_path}")
    
    def create_turbulence_properties(self):
        """Create turbulenceProperties file for laminar simulation"""
        print("\n=== Creating turbulenceProperties (laminar) ===")
        
        turbulence_props_content = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      turbulenceProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

simulationType laminar;

// ************************************************************************* //
"""
        
        turbulence_props_path = self.case_dir / "constant" / "turbulenceProperties"
        with open(turbulence_props_path, 'w') as f:
            f.write(turbulence_props_content)
        print(f"Created: {turbulence_props_path}")
    
    def create_initial_fields(self):
        """Create basic initial field files (U, T, p)"""
        print("\n=== Creating initial field files ===")
        
        # Velocity field (U)
        u_content = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 0 0);

boundaryField
{
    // Default boundary condition - modify as needed
    ".*"
    {
        type            zeroGradient;
    }
    
    // Example boundary conditions - uncomment and modify as needed
    /*
    inlet
    {
        type            fixedValue;
        value           uniform (100 0 0);  // Modify velocity as needed
    }
    
    outlet
    {
        type            zeroGradient;
    }
    
    ".*wall.*"
    {
        type            noSlip;
    }
    */
}

// ************************************************************************* //
"""
        
        u_path = self.case_dir / "0" / "U"
        with open(u_path, 'w') as f:
            f.write(u_content)
        print(f"Created: {u_path}")
        
        # Temperature field (T)
        t_content = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      T;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 1 0 0 0];

internalField   uniform 300;

boundaryField
{
    // Default boundary condition - modify as needed
    ".*"
    {
        type            zeroGradient;
    }
    
    // Example boundary conditions - uncomment and modify as needed
    /*
    inlet
    {
        type            fixedValue;
        value           uniform 600;  // Modify temperature as needed
    }
    
    outlet
    {
        type            zeroGradient;
    }
    
    ".*wall.*"
    {
        type            zeroGradient;  // or fixedValue for wall temperature
    }
    */
}

// ************************************************************************* //
"""
        
        t_path = self.case_dir / "0" / "T"
        with open(t_path, 'w') as f:
            f.write(t_content)
        print(f"Created: {t_path}")
        
        # Pressure field (p)
        p_content = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [1 -1 -2 0 0 0 0];

internalField   uniform 101325;

boundaryField
{
    // Default boundary condition - modify as needed
    ".*"
    {
        type            zeroGradient;
    }
    
    // Example boundary conditions - uncomment and modify as needed
    /*
    inlet
    {
        type            zeroGradient;  // or totalPressure
    }
    
    outlet
    {
        type            fixedValue;
        value           uniform 101325;  // Modify pressure as needed
    }
    
    ".*wall.*"
    {
        type            zeroGradient;
    }
    */
}

// ************************************************************************* //
"""
        
        p_path = self.case_dir / "0" / "p"
        with open(p_path, 'w') as f:
            f.write(p_content)
        print(f"Created: {p_path}")
    
    def run_preprocessing_commands(self):
        """Run necessary OpenFOAM preprocessing commands"""
        print("\n=== Running OpenFOAM preprocessing commands ===")
        
        commands = [
            ["checkMesh", "-case", str(self.case_dir)],
            ["changeDictionary", "-case", str(self.case_dir)]
        ]
        
        for cmd in commands:
            print(f"Running: {' '.join(cmd)}")
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                print(f"✓ {cmd[0]} completed successfully")
                if result.stdout:
                    print("Output:", result.stdout[-500:])  # Show last 500 chars
            except subprocess.CalledProcessError as e:
                print(f"⚠ Warning: {cmd[0]} failed: {e}")
                if e.stderr:
                    print("Error:", e.stderr[-500:])
    
    def run_complete_preprocessing(self):
        """Run the complete preprocessing workflow"""
        print("OpenFOAM NOZZLE PREPROCESSING WORKFLOW")
        print("=" * 50)
        
        try:
            # Step 1: Setup case directory
            self.setup_case_directory()
            
            # Step 2: Convert mesh
            self.convert_mesh_to_openfoam()
            
            # Step 3: Create system files
            self.create_control_dict()
            self.create_fv_schemes()
            self.create_fv_solution()
            self.create_change_dictionary_dict()
            
            # Step 4: Create constant files
            self.create_thermophysical_properties()
            self.create_turbulence_properties()
            
            # Step 5: Create initial fields
            self.create_initial_fields()
            
            # Step 6: Run preprocessing commands
            self.run_preprocessing_commands()
            
            print("\n" + "=" * 50)
            print("✓ OpenFOAM case setup completed successfully!")
            print(f"Case directory: {self.case_dir.absolute()}")
            print("\nNext steps:")
            print("1. Modify boundary conditions in 0/ files (U, T, p)")
            print("2. Adjust thermophysicalProperties if needed")
            print("3. Run changeDictionary if you modified changeDictionaryDict")
            print("4. Run: sonicFoam")
            print("=" * 50)
            
        except Exception as e:
            print(f"\n❌ Error during preprocessing: {e}")
            raise


def main():
    """Main function to run the OpenFOAM preprocessing"""
    # Initialize preprocessor
    preprocessor = OpenFOAMNozzlePreprocessor(
        mesh_file="structured_nozzle_3D.msh",
        case_name="nozzle_sonicfoam_case"
    )
    
    # Run complete preprocessing
    preprocessor.run_complete_preprocessing()


if __name__ == "__main__":
    main()
