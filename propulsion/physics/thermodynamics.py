from CoolProp.CoolProp import PropsSI

class n2oThermodynamics:
    def __init__(self, temperature):
        #parsed temperature is in Celsius
        self.temperature = temperature
        self.gas = 1
        self.liquid = 0
        self.fluidName = "NitrousOxide"

    def calculate_specific_gas_constant(self):
        return PropsSI('GAS_CONSTANT', self.fluidName) / PropsSI('M', self.fluidName)
    
    def all_properties(self, temperature=None):
        temp = temperature if temperature is not None else self.temperature
        density_vapour = PropsSI('D', 'T', temp, 'Q', self.gas, self.fluidName)  # promjena gustoce plina
        pressure_vapour = PropsSI('P', 'T', temp, 'Q', self.gas, self.fluidName)  # promjena tlaka plina

        density_liquid = PropsSI('D', 'T', temp, 'Q', self.liquid, self.fluidName)  # promjena gustoce tekucine
        pressure_liquid = PropsSI('P', 'T', temp, 'Q', self.liquid, self.fluidName)  # promjena tlaka tekucine
        return density_vapour, pressure_vapour, density_liquid, pressure_liquid
    
    def get_density_vapour(self, temperature=None):
        temp = temperature if temperature is not None else self.temperature
        return PropsSI('D', 'T', temp, 'Q', self.gas, self.fluidName)

    def get_pressure_vapour(self, temperature=None):
        temp = temperature if temperature is not None else self.temperature
        return PropsSI('P', 'T', self.temperature, 'Q', self.gas, self.fluidName)

    def get_density_liquid(self, temperature=None):
        temp = temperature if temperature is not None else self.temperature
        return PropsSI('D', 'T', self.temperature, 'Q', self.liquid, self.fluidName)

    def get_pressure_liquid(self, temperature=None):
        temp = temperature if temperature is not None else self.temperature
        return PropsSI('P', 'T', self.temperature, 'Q', self.liquid, self.fluidName)
    
    def get_temperature(self, pressure):
        return PropsSI('T', 'P', pressure, 'Q', self.liquid, self.fluidName)
    
    # def get_gamma_liquid(self, temperature=None):
    #     temp = temperature if temperature is not None else self.temperature
    #     #napraviti funkciju koja ce izracunati gamma
    #     Cp = PropsSI('C', 'T', temp, 'Q', self.liquid, self.fluidName) # J/kg/K
    #     print(Cp)
    #     R = PropsSI('GAS_CONSTANT', self.fluidName) / PropsSI('M', self.fluidName) #J/molK g/mol = J/gK
    #     R
    #     print(R)
    #     gamma = Cp/(Cp - R)
    #     return gamma

    def get_cp_liquid(self, temperature=None):
        temp = temperature if temperature is not None else self.temperature
        return PropsSI('C', 'T', temp, 'Q', self.liquid, self.fluidName)
    
    def get_enthalpy_liquid(self, temperature=None):
        temp = temperature if temperature is not None else self.temperature
        return PropsSI('H', 'T', temp, 'Q', self.liquid, self.fluidName)
    
    def get_enthalpy_vapour(self, temperature=None):
        temp = temperature if temperature is not None else self.temperature
        return PropsSI('H', 'T', temp, 'Q', self.gas, self.fluidName)

    def get_gamma_gas(self, pressure=101325, temperature=None): # napraviti tocniji izracun gamme
        # temp = temperature if temperature is not None else self.temperature
        # Cp = PropsSI('CPMOLAR', 'T', temp, 'Q', self.liquid, self.fluidName)
        # # Cv = PropsSI('CVMOLAR', 'T', temp, 'P', pressure, self.fluidName)
        # # return Cp / Cv
        # R = PropsSI('GAS_CONSTANT', self.fluidName) / PropsSI('M', self.fluidName) #J/molK g/mol = J/gK
        # # gamma = Cp/(Cp - R)
        # gamma = PropsSI('ISENTROPIC_EXPONENT', 'T', self.temperature, 'P', pressure, self.fluidName)
        gamma = 1.3
        return gamma




class fuel_thermodynamics:
    def __init__(self, temperature, pressure_liquid, FUEL='Acetone'):
        self.FUEL = FUEL
        self.pressure_liquid = pressure_liquid
        self.temperature = temperature

        # self.initial_fuel_density = self.get_density_fuel(self.temperature, self.pressure_liquid)

    # Functions
    def calculate_specific_gas_constant(self):
        return PropsSI('GAS_CONSTANT', self.FUEL) / PropsSI('M', self.FUEL)

    
    def get_density_fuel(self, pressure_liquid=None, temperature=None):
        if pressure_liquid is None:
            pressure_liquid = self.pressure_liquid
        if temperature is None:
            temperature = self.temperature
        return PropsSI('D', 'P', pressure_liquid, 'T', temperature, self.FUEL)

    def get_dynamic_viscosity(self, pressure_liquid=None, temperature=None):
        """
        Returns the dynamic viscosity of the fuel [Pa·s].
        """
        if pressure_liquid is None:
            pressure_liquid = self.pressure_liquid
        if temperature is None:
            temperature = self.temperature

        # print(pressure_liquid, temperature)
        return PropsSI('VISCOSITY', 'P', pressure_liquid, 'T', temperature, 'Ethanol')



# oxidizer = n2oThermodynamics(300)
# fuel = fuel_thermodynamics(300, oxidizer.get_pressure_liquid(), 'Acetone')
# print(fuel.get_density_fuel())