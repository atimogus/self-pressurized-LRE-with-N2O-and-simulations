from propulsion.engine.nozzleGeometry import NozzleGeometry
from propulsion.engine.scrintleInjector import Injector
from propulsion.physics.CEA import CEACantera
from propulsion.physics.thermodynamics import n2oThermodynamics, fuel_thermodynamics
from propulsion.tank.concentricTanks import ConcentricTanks
from simulations.HotFire import HotFireSimulation
from simulations.ColdFlowTest import ColdFlowTestSimulation
from CAD.regenCooling import regenNozzle
# from simulations.optimisationCFD import optimisationCFD
from simulations.CFD.ngo.nozzleGeometryOptimisation import start_optimisation
import json
import os
import subprocess

from config import Config

class initial_design:

    def __init__(self):
        self.oxidizer = n2oThermodynamics(Config.INITIAL_TEMP_K)
        self.fuel = fuel_thermodynamics(Config.INITIAL_TEMP_K, self.oxidizer.get_pressure_liquid(), Config.FUEL)

        self.tank = ConcentricTanks(
            oxidizer_mass_kg=Config.OXIDIZER_MASS_KG,
            oxidizer_density_liquid=self.oxidizer.get_density_liquid(),
            oxidizer_density_vapour=self.oxidizer.get_density_vapour(),
            density_fuel=self.fuel.get_density_fuel(),
            o_f_ratio=Config.OXIDIZER_TO_FUEL_RATIO,
            fuel_duration=Config.FUEL_DURATION,
            tank_material_strength=Config.TANK_MATERIAL_STRENGTH,
            safety_factor=Config.SAFETY_FACTOR_TANK,
            tank_ullage = Config.TANK_ULAGE
        )
        self.tank.print_tank_parameters()

        self.scrintleInjector = Injector(
            discharge_coefficient=Config.DISCHARGE_COEFFICIENT,
            discharge_coefficient_fuel=Config.discharge_coefficient_fuel,
            number_of_fuel_orfices=Config.NUMBER_OF_ORIFICES_FUEL,
            numberOf_oxidizer_orfices=Config.NUMBER_OF_ORIFICES_OXIDIZER,
            OXIDIZER_ORIFICE_DIAMETER=Config.OXIDIZER_ORIFICE_DIAMETER,
            oxidizer_liquid_density=self.oxidizer.get_density_liquid(),
            oxidizer_liquid_pressure=self.oxidizer.get_pressure_liquid(),
            combustion_pressure=Config.COMBUSTION_CHAMBER_PRESSURE,
            oxidizer_fuel_ratio=Config.OXIDIZER_TO_FUEL_RATIO,
            temperature=Config.INITIAL_TEMP_K,
            density_f=self.fuel.get_density_fuel(),
            PRESSURE_FRICTION_LOSS=Config.PRESSURE_FRICTION_LOSS)


        # Initialize CEACantera object
        self.cea = CEACantera(Config.YAML_FILE, Config.FUEL_FOR_CEA, Config.OXIDIZER, Config.INITIAL_TEMP_K, 
                                Config.COMBUSTION_CHAMBER_PRESSURE, Config.OXIDIZER_TO_FUEL_RATIO)
        self.cea.resoultsCEA()

        # Initialize NozzleGeometry object
        self.nozzle = NozzleGeometry(config = Config,
                                    cea=self.cea,
                                    injector = self.scrintleInjector)
        self.nozzle.printNozzleGeometry()

        self.regenNozzleDesign = regenNozzle(
            combustion_chamber_diameter=self.nozzle.calculate_combustion_chamber_diameter() + 2 * Config.CHAMBERSAFE_THICKNESS/1000,
            vanjski_promjer_ulaz=Config.OXIDIZER_TANK_DIAMETER,
            promjer_grla_mlaznice=self.nozzle.circle_diameter(self.nozzle.calculate_throat_area()) * 1000,
            tlak_komore=Config.COMBUSTION_CHAMBER_PRESSURE,
            tlak_fluida=self.oxidizer.get_pressure_liquid(),
            cvrstoca_materijala=Config.TANK_MATERIAL_STRENGTH,
            faktor_sigurnosti=Config.SAFETY_FACTOR_CHAMBER,
            d_hidraulicke_cijevi=Config.D_HIDRAULIC_HOSE,
            maseni_protok_IPA=self.scrintleInjector.fuel_massflowrate,
            gustoca_IPA=785,
            alu_thickness=Config.ALU_THICKNESS,
            n=Config.NUNMER_OF_REGEN_CHANNELS,
            k=1
        )

        self.regenNozzleDesign.printDesign()

        # self.optimisation_CFD = optimisationCFD(config = Config,
        #             fuel = self.fuel, 
        #             injector = self.scrintleInjector,
        #             nozzle = self.nozzle,
        #             oxidizer = self.oxidizer
        #             )

        # self.optimisation_CFD.print_results()

        if Config.execute_simulation:
            if Config.Cold_Flow_test:
                # Initialize ColdFlowTestSimulation object
                self.simulation = ColdFlowTestSimulation(
                    oxidizer = self.oxidizer,
                    fuel = self.fuel,
                    config=Config,
                    cea=self.cea,
                    tank=self.tank,
                    injector=self.scrintleInjector,
                    nozzle=self.nozzle,
                )
            else:
                # Initialize HotFireSimulation object
                self.simulation = HotFireSimulation(
                    oxidizer = self.oxidizer,
                    fuel = self.fuel,
                    config=Config,
                    cea=self.cea,
                    tank=self.tank,
                    injector=self.scrintleInjector,
                    nozzle=self.nozzle,
                )



if __name__ == "__main__":
    initialMotorDesign = initial_design()

    if Config.execute_simulation:
        initialMotorDesign.simulation.run_simulation()
        initialMotorDesign.simulation.plot_simulation_results()

        if Config.flight_simulation and not Config.Cold_Flow_test:
            print("Flight simulation is not implemented yet.")


    if Config.generate_Freecad_geometry:
        # Step 1: Extract data from regenNozzleDesign
        regen_data = {
            "fcstd_filename": "mlaznicaCFD.FCStd",

            "D_inlet": initialMotorDesign.regenNozzleDesign.D_inlet,
            "D_ulaz_mlaznice": initialMotorDesign.regenNozzleDesign.D_cooling_channels_inlet,
            "d_cooling_channels_inlet": initialMotorDesign.regenNozzleDesign.d_cooling_channels_inlet,
            "combustion_chamber_diameter": initialMotorDesign.regenNozzleDesign.combustion_chamber_diameter,

            "promjer_grla_mlaznice": initialMotorDesign.regenNozzleDesign.promjer_grla_mlaznice,
            "d_cooling_channels_throat": initialMotorDesign.regenNozzleDesign.d_cooling_channels_throat,
            "D_cooling_channels_throat": initialMotorDesign.regenNozzleDesign.D_cooling_channels_throat,
            "D_throat": initialMotorDesign.regenNozzleDesign.D_throat,

            "fi_revolve": initialMotorDesign.regenNozzleDesign.fi_revolve,
            "fi_alu_stupnjevi": initialMotorDesign.regenNozzleDesign.fi_alu_stupnjevi,
            "fi_inlet_stupnjevi_ulaz": initialMotorDesign.regenNozzleDesign.fi_inlet_stupnjevi_ulaz,
            "fi_inlet_stupnjevi_izlaz": initialMotorDesign.regenNozzleDesign.fi_inlet_stupnjevi_izlaz,
            "channel_height": initialMotorDesign.regenNozzleDesign.channel_height,

            "CONVERGENT_NOZZLE_ANGLE": Config.CONVERGENT_NOZZLE_ANGLE,
            "DIVERGENT_NOZZLE_ANGLE": Config.DIVERGENT_NOZZLE_ANGLE,
            "Combustion_Chamber_Height": initialMotorDesign.nozzle.combustion_chamber_height,
            "Divergent_Nozzle_Length": initialMotorDesign.nozzle.divergent_nozzle_length,
            "Convergent_Nozzle_Length": initialMotorDesign.nozzle.convergent_nozzle_length,
        }

        # Step 2: Save the data to a JSON file
        data_file = os.path.join(os.path.dirname(__file__), "CAD", "design_data.json")
        with open(data_file, "w") as f:
            json.dump(regen_data, f)

        print(f"\nDesign data saved to {data_file}")

        # Step 3: Start FreeCAD with the geometryFreecad.py script
        geometry_script = os.path.join(os.path.dirname(__file__), "CAD", "geometryFreecad.py")
        freecad_appimage = os.path.abspath("./CAD/FreeCAD1.0.AppImage")  # Full path to the FreeCAD AppImage
        subprocess.run([freecad_appimage, geometry_script])

    if Config.nozzle_contour_optimization:
        start_optimisation(y_combustion_chamber = initialMotorDesign.nozzle.combustion_chamber_diameter / 2,
                           y_exit = initialMotorDesign.nozzle.exit_diameter / 2,
                           y_throat = initialMotorDesign.nozzle.throat_diameter / 2,
                           num_of_points_nozzle = Config.num_of_points_nozzle
                           )
        
                           

