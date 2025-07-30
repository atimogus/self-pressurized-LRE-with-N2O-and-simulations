import numpy as np

#todo: calculate calculate_combustion_chamber_height based on injector settings (time_spray+ time_combustion + time_mixing)

class NozzleGeometry:
    def __init__(self, config, cea, injector):
        self.config = config
        self.cea = cea
        self.injector = injector

        self.combinedMassFlowRate = injector.oxidizer_massflowrate + injector.fuel_massflowrate

        # Calculate nozzle parameters
        self.exit_velocity = self.calculate_exit_velocity()
        self.throat_area = self.calculate_throat_area()
        self.exit_area = self.calculate_exit_area()
        self.aspect_ratio = self.calculate_aspect_ratio()
        self.isp = self.calculate_isp()
        self.thrust = self.calculate_thrust()
        self.throat_diameter = self.circle_diameter(self.throat_area)
        self.exit_diameter = self.circle_diameter(self.exit_area)
        self.combustion_chamber_diameter = self.calculate_combustion_chamber_diameter()
        self.combustion_chamber_height = self.calculate_combustion_chamber_height()
        self.divergent_nozzle_length = self.calculate_divergent_nozzle_length()
        self.convergent_nozzle_length = self.calculate_convergent_nozzle_length()

    def calculate_exit_velocity(self):
        Ve = np.sqrt(
            2 * self.cea.gam * self.cea.R * self.cea.gas.T / (self.cea.gam - 1) *
            (1 - (self.config.P_EXIT / self.config.COMBUSTION_CHAMBER_PRESSURE)**((self.cea.gam - 1) / self.cea.gam))
        )
        return Ve

    def calculate_throat_area(self):
        At =  self.combinedMassFlowRate / ((self.config.COMBUSTION_CHAMBER_PRESSURE * np.sqrt(self.cea.gam)) /
                                            np.sqrt(self.cea.R * self.cea.gas.T) *
                                              (2 / (self.cea.gam + 1))**((self.cea.gam + 1) / (2 * (self.cea.gam - 1))))
        return At

    def calculate_exit_area(self):
        At = self.calculate_throat_area()
        # Ae = At * (
        #     (self.config.P_EXIT / self.config.COMBUSTION_CHAMBER_PRESSURE)**(-1 / self.cea.gam) *
        #     np.sqrt((1 - (self.config.P_EXIT / self.config.COMBUSTION_CHAMBER_PRESSURE)**((self.cea.gam - 1) / self.cea.gam)) / ((self.cea.gam - 1) / 2))
        # )
        Ae = At * np.sqrt((self.cea.gam - 1) / 2) * (2 / (self.cea.gam + 1))**((self.cea.gam + 1) / (2 * (self.cea.gam - 1))) * 1 / ((self.config.P_EXIT / self.config.COMBUSTION_CHAMBER_PRESSURE)**(1 / self.cea.gam) * np.sqrt(1 - (self.config.P_EXIT / self.config.COMBUSTION_CHAMBER_PRESSURE)**((self.cea.gam - 1) / self.cea.gam)))
        return Ae

    def circle_diameter(self, area):
        diameter = np.sqrt(4 * area / np.pi)
        return diameter

    def calculate_combustion_chamber_diameter(self):
        throat_diameter = self.circle_diameter(self.calculate_throat_area())
        combustion_chamber_diameter = 3.5 * throat_diameter
        return combustion_chamber_diameter

    def calculate_aspect_ratio(self):
        At = self.calculate_throat_area()
        Ae = self.calculate_exit_area()
        aspectRatio = Ae / At
        return aspectRatio

    def calculate_isp(self):
        Ve = self.calculate_exit_velocity()
        Isp = Ve / self.config.gravity_constant
        return Isp

    def calculate_thrust(self):
        Ve = self.calculate_exit_velocity()
        Ae = self.calculate_exit_area()
        seaLevelThrust = self.combinedMassFlowRate * Ve + Ae * (self.config.P_EXIT - self.config.P_ATMOSPHERE)
        return seaLevelThrust

    def calculate_combustion_chamber_height(self):
        # Assuming a cylindrical combustion chamber, it must be calculated on injector settings
        combustion_chamber_diameter = self.calculate_combustion_chamber_diameter()
        combustion_chamber_height = 1.2 * combustion_chamber_diameter
        return combustion_chamber_height

    def calculate_divergent_nozzle_length(self):
        # Assuming a simple conical nozzle
        divergent_nozzle_length = (self.exit_diameter/2 - self.throat_diameter/2) / np.tan(np.radians(16.6))
        return divergent_nozzle_length

    def calculate_convergent_nozzle_length(self):
        # Assuming a simple conical nozzle
        convergent_nozzle_length = (self.combustion_chamber_diameter/2 - self.throat_diameter/2 ) / np.tan(np.radians(45))
        return convergent_nozzle_length
    
    #@staticmethod
    def printNozzleGeometry(self):
        print("\n---------------------------- INITIAL Nozzle Geometry:-------------------------------")
        # Print nozzle parameters with values rounded to the second decimal
        print(f"Nozzle Exit Velocity: {self.exit_velocity:.2f} m/s")
        print(f"Throat Area: {self.throat_area} m^2")
        print(f"Exit Area: {self.exit_area} m^2")
        print(f"Aspect Ratio: {self.aspect_ratio:.2f}")
        print(f"Specific Impulse: {self.isp:.2f} s")
        print(f"Thrust: {self.thrust:.2f} N")
        print(f"Throat Diameter: {self.throat_diameter * 1000:.2f} mm")
        print(f"Exit Diameter: {self.exit_diameter * 1000:.2f} mm")
        print(f"Combustion Chamber Diameter: {self.combustion_chamber_diameter * 1000:.2f} mm")
        print(f"Combustion Chamber Height: {self.combustion_chamber_height * 1000:.2f} mm")
        print(f"Divergent Nozzle Length: {self.divergent_nozzle_length * 1000:.2f} mm")
        print(f"Convergent Nozzle Length: {self.convergent_nozzle_length * 1000:.2f} mm")