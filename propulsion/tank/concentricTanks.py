import numpy as np

# todo: edit this class so comfig file is passed and fuel and oxidizer, not this values
class ConcentricTanks:
    def __init__(self, oxidizer_mass_kg: float, oxidizer_density_liquid: float, 
                 oxidizer_density_vapour: float, density_fuel: float, 
                 o_f_ratio: float, fuel_duration: float, 
                 tank_material_strength: float, safety_factor: float, tank_ullage):
        
        self.oxidizer_mass_kg = oxidizer_mass_kg
        self.oxidizer_density_liquid = oxidizer_density_liquid
        self.oxidizer_density_vapour = oxidizer_density_vapour
        self.density_fuel = density_fuel
        self.o_f_ratio = o_f_ratio
        self.fuel_duration = fuel_duration
        self.tank_material_strength = tank_material_strength
        self.safety_factor = safety_factor
        self.tank_ullage = tank_ullage  # Default value for ullage (can be parameterized)

        self.mass_fuel = self.initial_fuel_mass()
        self.volume_tank_oxidizer = self.initial_oxidizer_volume()[2]
        self.mass_vapour_ox = self.oxidizer_initial_mass_vapour()


    def initial_oxidizer_volume(self):
        volume_liquid = self.oxidizer_mass_kg / self.oxidizer_density_liquid
        volume_gas = self.tank_ullage * volume_liquid
        volume_tank_oxidizer = volume_liquid + volume_gas
        return volume_liquid, volume_gas, volume_tank_oxidizer

    def oxidizer_initial_mass_vapour(self):
        _, volume_gas, _ = self.initial_oxidizer_volume()
        oxidizer_initial_mass_vapour = volume_gas * self.oxidizer_density_vapour
        return oxidizer_initial_mass_vapour

    def initial_fuel_mass(self):
        mass_fuel = (self.oxidizer_mass_kg * self.fuel_duration)/ self.o_f_ratio 
        return mass_fuel

    def initial_fuel_volume(self):
        fuel_volume = self.mass_fuel / self.density_fuel
        return fuel_volume

    def tank_dimensions_l_d(self, l_d_ratio):
        # Calculating dimensions with l/d ratio
        volume_fuel = self.initial_fuel_volume()
        volume_liquid, volume_gas, _ = self.initial_oxidizer_volume()
        volume_tank = volume_fuel + volume_gas + volume_liquid

        body_diameter = ((volume_tank * 4) / (np.pi * l_d_ratio)) ** (1 / 3)
        tank_height = l_d_ratio * body_diameter
        fuel_tank_diameter = np.sqrt((volume_fuel * 4) / (tank_height * np.pi))

        rocket_body_area = (body_diameter ** 2) * np.pi / 4
        return body_diameter, tank_height, fuel_tank_diameter, rocket_body_area

    def tank_dimensions_via_diameter(self, F_TANK_DIAM, OX_TANK_DIAM, OX_TANK_THICKNESS, F_TANK_THICKNESS):
        # Calculate dimensions if Diameters of both tanks are given
        volume_fuel = self.initial_fuel_volume()
        volume_liquid, volume_gas, _ = self.initial_oxidizer_volume()
        
        # body_diameter = ((Volume_tank * 4) / (3.14 * L_D_ratio))**(1/3)
        # tank_height = L_D_ratio * body_diameter
        # fuel_tank_diameter = np.sqrt((V_fuel*4)/(tank_height * 3.14))
        
        # ROCKET_BODY_AREA = (body_diameter) ** 2 * np.pi / 4

        f_tank_height_via_dia = (4 * volume_fuel )/(np.pi * ((F_TANK_DIAM - 2 * F_TANK_THICKNESS) / 1000)**2)
        ox_tank_height_via_dia = (4 * (volume_liquid+ volume_gas))/(np.pi * (((OX_TANK_DIAM - 2 * OX_TANK_THICKNESS) / 1000)**2 - ((F_TANK_DIAM) / 1000)**2))
        ox_tank_height_via_dia_without_gas = (4 * volume_liquid)/(np.pi * (((OX_TANK_DIAM - 2 * OX_TANK_THICKNESS) / 1000)**2 - ((F_TANK_DIAM) / 1000)**2))

        return f_tank_height_via_dia, ox_tank_height_via_dia , ox_tank_height_via_dia_without_gas


    def tank_thickness(self):
        body_diameter, _, _, _ = self.tank_dimensions()
        tmin1 = (self.oxidizer_density_vapour * body_diameter / 2) / (self.tank_material_strength / self.safety_factor)
        tmin2 = (self.oxidizer_density_vapour * body_diameter / 2) / (2 * self.tank_material_strength / self.safety_factor)

        if tmin1 >tmin2:
            return tmin1
        else:
            return tmin2
        
    def print_tank_parameters(self):
        print(f"\n------------------------Tank Parameters=============================")
        print(f"self.mass_fuel: {self.mass_fuel}")
        print(f"oxidizer_mass_kg: {self.oxidizer_mass_kg}")
        print(f"self.mass_vapour_ox: {self.mass_vapour_ox}")



    # def gas_orfice_sizing(fuel_tank_diameter, tank_height):
    #     # Calculate the orifice diameter in fuel tank so vapour can pass from oxigen tank to fuel tank
    #     sim_time = np.array(data["pressure_c"]) / TIME_STEP
    #     sim_time_sequence = np.arange(0, len(sim_time) * TIME_STEP, TIME_STEP)
    #     # print(sim_time_sequence[-1], fuel_tank_diameter, tank_height)
    #     return np.sqrt((fuel_tank_diameter**2 * tank_height)/(ORFICE_NUMBER * GAS_VELOCIY * sim_time_sequence[-1]))