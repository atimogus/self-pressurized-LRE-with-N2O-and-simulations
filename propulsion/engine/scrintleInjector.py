import numpy as np

# TODO:
# chemical kinetics, how to determine chamber colume
# HEM and DYER methods for mass flow rate

class Injector:
    """Class to model the behavior of a rocket engine injector."""

    def __init__(self, discharge_coefficient, discharge_coefficient_fuel, number_of_fuel_orfices, numberOf_oxidizer_orfices, OXIDIZER_ORIFICE_DIAMETER, oxidizer_liquid_density, 
                 oxidizer_liquid_pressure, combustion_pressure, oxidizer_fuel_ratio, temperature,
                 density_f, PRESSURE_FRICTION_LOSS):
        self.discharge_coefficient = discharge_coefficient # discharge coefficient is the same for fuel and oxidizer (in reality this is not the case!!!!)
        self.discharge_coefficient_fuel = discharge_coefficient_fuel
        self.number_of_fuel_orfices = number_of_fuel_orfices
        self.numberOf_oxidizer_orfices = numberOf_oxidizer_orfices
        self.OXIDIZER_ORIFICE_DIAMETER = OXIDIZER_ORIFICE_DIAMETER / 1000 # convert mm to m
        self.oxidizer_liquid_density = oxidizer_liquid_density
        self.oxidizer_liquid_pressure = oxidizer_liquid_pressure
        self.combustion_pressure = combustion_pressure
        self.oxidizer_fuel_ratio = oxidizer_fuel_ratio
        self.temperature = temperature
        self.density_f = density_f
        self.PRESSURE_FRICTION_LOSS = PRESSURE_FRICTION_LOSS

        self.areaOf_oxidizer_orfice = self.calculate_areaOf_oxidizer_orfice(self.OXIDIZER_ORIFICE_DIAMETER)
        self.oxidizer_massflowrate = self.calculate_mass_flow_rate_SPI(self.oxidizer_liquid_density, self.oxidizer_liquid_pressure, self.combustion_pressure)
        self.fuel_massflowrate = self.calculate_fuel_massflowrate(self.oxidizer_massflowrate)
        self.area_ofFuel_orfice = self.calculate_fuel_orifice_diameter()[0]

    def calculate_areaOf_oxidizer_orfice(self, diameter):
        # Calculate the area of the oxidizer orifice
        # A = (pi/4) * D^2
        areaOf_oxidizer_orfice = (np.pi / 4) * diameter**2
        return areaOf_oxidizer_orfice

    def calculate_mass_flow_rate_SPI(self, oxidizer_liquid_density, oxidizer_liquid_pressure, combustion_pressure):
        # assumption - single phase flow
        # Calculate the mass flow rate of the oxidizer and fuel using the SPI method
        # https://ntrs.nasa.gov/api/citations/20190001326/downloads/20190001326.pdf

        oxidizer_massflowRate = self.numberOf_oxidizer_orfices * self.discharge_coefficient * self.areaOf_oxidizer_orfice * np.sqrt(
            2 * oxidizer_liquid_density * (oxidizer_liquid_pressure - combustion_pressure))
        
        return oxidizer_massflowRate

    def calculate_fuel_massflowrate(self, oxidizer_massflowrate): 
        return oxidizer_massflowrate / self.oxidizer_fuel_ratio 

    def calculate_fuel_massflowrate_SPI(self, density_f, fuel_pressure, combustion_pressure):
        #dp is pressure delta = pressure higher - pressure
        fuel_massflowRate = self.number_of_fuel_orfices * self.discharge_coefficient_fuel * self.area_ofFuel_orfice * np.sqrt(2 * density_f * (fuel_pressure - combustion_pressure))
        
        
        return fuel_massflowRate

    def calculate_mass_flow_rate_HEM(self):
        # assumption - two phase flow
        # https://ntrs.nasa.gov/api/citations/20190001326/downloads/20190001326.pdf
        return 0 

    def calculate_mass_flow_rate_DYER (self):
        # assumption - two phase flow - better approximation than HEM
        # https://ntrs.nasa.gov/api/citations/20190001326/downloads/20190001326.pdf
        return 0 

    def calculate_fuel_orifice_diameter(self):
        A_fuel = self.fuel_massflowrate / (
            self.number_of_fuel_orfices * self.discharge_coefficient_fuel * np.sqrt(2 * self.density_f * (self.oxidizer_liquid_pressure - self.PRESSURE_FRICTION_LOSS - self.combustion_pressure))
        )
        Diameter_fuel_m = np.sqrt(A_fuel * 4 / np.pi)
        return A_fuel, Diameter_fuel_m

    def calculate_ejection_velocity(self):
        A_fuel, _ = self.calculate_fuel_orifice_diameter()
        speed_OxInjector = self.oxidizer_massflowrate / (self.oxidizer_liquid_density * self.areaOf_oxidizer_orfice)
        speed_FInjector = self.fuel_massflowrate / (self.density_f * A_fuel)
        return speed_OxInjector, speed_FInjector