import cantera as ct

class CEACantera:
    def __init__(self, yaml_file, fuel, oxidizer, init_temp, pressure_chamber, OXIDIZER_TO_FUEL_RATIO):
        """
        Initialize the CEACantera object with the required parameters.
        """
        self.yaml_file = yaml_file
        self.fuel = fuel
        self.oxidizer = oxidizer
        self.init_temp = init_temp
        self.pressure_chamber = pressure_chamber
        self.gas = ct.Solution(yaml_file)
        self.OXIDIZER_TO_FUEL_RATIO = OXIDIZER_TO_FUEL_RATIO

        self.calculate(self.OXIDIZER_TO_FUEL_RATIO)

        self.prandtl_number = self.gas.viscosity * self.gas.cp_mass / self.gas.thermal_conductivity
    

    def calculate(self, o_f, temperature=None, pressure_chamber=None):
        """
        Perform equilibrium calculations for the given oxidizer-to-fuel ratio (o_f).
        """
        # Use default values if temperature or pressure_chamber are not provided
        if temperature is None:
            temperature = self.init_temp
        if pressure_chamber is None:
            pressure_chamber = self.pressure_chamber

        # Calculate mass fractions
        m_fuel = 1 / (1 + o_f)  # Fuel mass fraction
        m_ox = o_f / (1 + o_f)  # Oxidizer mass fraction

        # Set the state of the gas mixture
        self.gas.TPY = temperature, pressure_chamber, {self.fuel: m_fuel, self.oxidizer: m_ox}
        
        # Perform equilibrium calculation at constant pressure and enthalpy
        self.gas.equilibrate('HP')
        
        # Calculate additional properties
        self.gam = self.gas.cp_mass / self.gas.cv_mass  # Specific heat ratio
        self.R = ct.gas_constant / self.gas.mean_molecular_weight  # Gas constant
    
    def resoultsCEA(self):
        print("\n---------------------------- CEA Results:-------------------------------")
        print(f"Temperature: {self.gas.T:.2f} K")
        print(f"Density: {self.gas.density:.2f} kg/m^3")
        print(f"Enthalpy: {self.gas.enthalpy_mass:.2f} J/kg")
        print(f"Gibbs: {self.gas.gibbs_mass:.2f} J/kg")
        print(f"Entropy: {self.gas.entropy_mass:.2f} J/kg")
        print(f"Mean Molecular Weight: {self.gas.mean_molecular_weight:.2f} g/mol")
        print(f"Specific Heat at Constant Pressure: {self.gas.cp_mass:.2f} J/(kg*K)")
        print(f"Specific Heat Ratio: {self.gam:.2f}")
        print(f"Gas Constant: {self.R:.2f} J/(kg*K)")
        print(f"Viscosity: {self.gas.viscosity} Pa*s")
        print(f"Thermal Conductivity: {self.gas.thermal_conductivity} W/(m*K)")
        print(f"Prandtl Number: {self.prandtl_number:.2f}")