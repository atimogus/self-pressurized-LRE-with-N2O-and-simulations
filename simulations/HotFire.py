import numpy as np
import matplotlib.pyplot as plt

'''napraviti evaporaciju prema modelu iz rada : https://arxiv.org/pdf/2302.06725
    stranica 21 - 22
'''

class HotFireSimulation:
    def __init__(self, oxidizer, fuel, config, cea, tank, injector, nozzle):
        self.config = config
        self.oxidizer = oxidizer
        self.cea = cea
        self.fuel = fuel
        self.tank = tank
        self.injector = injector
        self.nozzle = nozzle

        self.ox_pressure_liquid = [] # pressure liquid for oxidizer is the same as pressure of vapour in tank beacause it is looked like it is in every timestep in equilibrium (saturation pressure)
        self.ox_liquid_massflowrate = []
        self.ox_vapour_massflowrate = []
        self.mass_liquid_ox = []
        self.mass_vapour_ox = []
        self.mass_fuel = []
        self.fuel_massflowrate = []
        self.O_F_shift = []
        self.fuel_pressure = []
        self.dp_dt = []             # this value is used to calculate the pressure drop in the tank, in real world this can be approximated as constant (done in this code by taking the mean of the last 5 values)
        self.pressure_combustion = []
        self.temperature_combustion = []
        self.total_thrust = []
        self.density_exhaust_gas = []
        self.temperature_tank = []
        self.rockett_acceleration = []
        self.rocket_velocity = []
        self.rocket_distance = []
        self.TW_ratio = []
        self.time = []

        # Initialize state variables
        self.ox_pressure_liquid.append(self.injector.oxidizer_liquid_pressure)
        self.temperature_tank.append(self.config.INITIAL_TEMP_K)
        self.mass_liquid_ox.append(self.tank.oxidizer_mass_kg)
        self.ox_liquid_massflowrate.append(self.injector.oxidizer_massflowrate * self.config.TIME_STEP)

        self.mass_fuel.append(self.tank.mass_fuel)
        self.fuel_massflowrate.append(self.injector.fuel_massflowrate * self.config.TIME_STEP)
        self.mass_vapour_ox.append(self.tank.mass_vapour_ox)

        self.pressure_combustion.append(self.config.COMBUSTION_CHAMBER_PRESSURE)
        self.fuel_pressure.append(self.injector.oxidizer_liquid_pressure - self.config.PRESSURE_FRICTION_LOSS)
        self.total_thrust.append(self.nozzle.thrust)


    def run_simulation(self):

        iteration = 0
        pressure_old = 0
        o_f_old = 0

        while self.mass_liquid_ox[iteration] > self.ox_liquid_massflowrate[iteration]:

            # Update the simulation state
            pressure_old, o_f_old = self.run_liquid_phase_simulation(iteration, pressure_old, o_f_old)
            self.update_time(iteration)
            iteration += 1

        mass_vapour_ox_initial = self.mass_vapour_ox[iteration]

        while self.mass_vapour_ox[iteration] > mass_vapour_ox_initial * 0.05 and self.ox_pressure_liquid[iteration] > 2* self.config.P_ATMOSPHERE: # pressure liquid is the same as pressure of vapour because it was all time in equilibrium (saturation pressure)
            self.run_vapor_phase_simulation(iteration)
            self.update_time(iteration)
            iteration += 1
        
        self.update_time(iteration)


    def run_liquid_phase_simulation(self, iteration, pressure_old, o_f_old):
        #writing temperature this way so it can be overwritten later in code
        temperature = self.temperature_tank[iteration]

        # Calculate the properties of the oxidizer
        ox_density_vapour, _, ox_density_liquid, ox_pressure_liquid = self.oxidizer.all_properties( temperature )

        #this section is needed for simulation stability, if finetune_value is set to 0 you can see instabilities in simulation
        finetune_value = 1
        if self.mass_liquid_ox[iteration] < self.config.OXIDIZER_MASS_KG * 0.5 and (
            pressure_old - ox_pressure_liquid < np.mean(self.dp_dt[-finetune_value:]) * 0.5):

            if self.mass_fuel[iteration] > self.fuel_massflowrate[iteration]:
                finetune_value = 1 #smaller numebr is better because it is less of constraint on the simulation
                temperature = self.oxidizer.get_temperature(pressure_old - np.mean(self.dp_dt[-finetune_value:]))
                ox_density_vapour, _, ox_density_liquid, ox_pressure_liquid = self.oxidizer.all_properties( temperature)

            else:
                finetune_value = 1  #smaller numebr is better because it is less of constraint on the simulation
                temperature = self.oxidizer.get_temperature(pressure_old - np.mean(self.dp_dt[-finetune_value:]))
                ox_density_vapour, _, ox_density_liquid, ox_pressure_liquid = self.oxidizer.all_properties( temperature)

        #calculating massflowrate of oxidizer
        ox_liquid_massflowrate = self.injector.calculate_mass_flow_rate_SPI(ox_density_liquid, ox_pressure_liquid, self.pressure_combustion[iteration]) * self.config.TIME_STEP
        mass_liquid_ox_new_it = self.mass_liquid_ox[iteration] - ox_liquid_massflowrate 
        mass_liquid_ox_old = mass_liquid_ox_new_it
        # print(" ox_liquid_massflowrate ", ox_liquid_massflowrate, " ox_density_liquid ",ox_density_liquid, " ox_pressure_liquid ", ox_pressure_liquid, " self.pressure_combustion[iteration] ", self.pressure_combustion[iteration] )
        # self.tank.volume_tank_oxidizer += massflowrate_fuel / density_fuel #povecanje volumena koj plin mora zauzeti
        # mass_liquid_ox_new_it, mass_vapour_ox, mass_vapourised = self.mass_evapouration(
        #     self.tank.volume_tank_oxidizer, mass_liquid_ox_old, self.mass_vapour_ox[iteration], ox_density_vapour, ox_density_liquid)


        if self.mass_fuel[iteration] > 0:

            # print("mass_fuel", self.mass_fuel[iteration], round(self.mass_fuel[iteration], 2), iteration)

            #calculating massflowrate of oxidizer with pressure drop because of friction in piston
            pressure_fuel = ox_pressure_liquid - self.config.PRESSURE_FRICTION_LOSS
            density_fuel = self.fuel.get_density_fuel(pressure_fuel, temperature)
            massflowrate_fuel = self.injector.calculate_fuel_massflowrate_SPI( density_fuel, pressure_fuel, self.pressure_combustion[iteration] ) * self.config.TIME_STEP
            # print("massflowrate_fuel", massflowrate_fuel, "pressure_fuel", pressure_fuel, "density_fuel", density_fuel, "self.pressure_combustion[iteration]", self.pressure_combustion[iteration])
            mass_fuel = self.mass_fuel[iteration] - massflowrate_fuel # getting new mass of fuel in tank after mass of fuel is gone through injector

            self.tank.volume_tank_oxidizer += massflowrate_fuel / density_fuel #povecanje volumena koj plin mora zauzeti
            mass_liquid_ox_new_it, mass_vapour_ox, mass_vapourised = self.mass_evapouration(
                self.tank.volume_tank_oxidizer, mass_liquid_ox_old, self.mass_vapour_ox[iteration], ox_density_vapour, ox_density_liquid)

            o_f = ox_liquid_massflowrate / massflowrate_fuel
            # print("o_f", o_f, "massflowrate_fuel", massflowrate_fuel, "ox_liquid_massflowrate", ox_liquid_massflowrate)

            if np.isnan(o_f):
                print(f"\nERROR: o_f is NaN, massflowrate_fuel: {massflowrate_fuel}, ox_liquid_massflowrate: {ox_liquid_massflowrate}")
                print(f"lowest pressure in tank: {pressure_fuel}, pressure_combustion: {self.pressure_combustion[iteration]}")  
                print(f"temperature: {temperature}\n")
                raise ValueError("ERROR: your chamber pressure is too high or your initial temperature is too low")
            else:
                if int(o_f * 100) % 10 != int(o_f_old * 100) % 10:
                    self.cea.calculate(o_f, temperature, self.pressure_combustion[iteration])
                

            thrust = (ox_liquid_massflowrate + massflowrate_fuel) / self.config.TIME_STEP * self.calculate_escape_velocity(
                self.cea.gam, self.cea.gas.T, self.cea.gas.mean_molecular_weight, self.config.P_ATMOSPHERE, self.pressure_combustion[iteration]) # nije P_EXIT ni P_ATMOSPHERE zapravo jer se tlak mjenja u kako raketa leti, potrebno integirrati flight dynamics za bolju tocnost
            
            pressure_combustion_new = self.calculate_new_combustion_pressure(
                        self.cea.gam, self.cea.gas.T, (ox_liquid_massflowrate + massflowrate_fuel) / self.config.TIME_STEP, 
                        self.nozzle.throat_area, self.cea.R
                    )
            o_f_old = o_f
        
            self.fuel_pressure.append(pressure_fuel)
            self.O_F_shift.append(o_f)
            self.fuel_massflowrate.append(massflowrate_fuel)
            self.pressure_combustion.append(pressure_combustion_new)
            self.density_exhaust_gas.append(self.cea.gas.density)

        else:
            self.pressure_combustion.append(self.config.P_ATMOSPHERE) # nije P_EXIT ni P_ATMOSPHERE zapravo jer se tlak mjenja u kako raketa leti, potrebno integirrati flight dynamics za bolju tocnost
            massflowrate_fuel = 0
            mass_fuel = 0
            thrust = 0
            self.fuel_massflowrate.append(massflowrate_fuel)
            self.fuel_pressure.append(0)

            mass_liquid_ox_new_it, mass_vapour_ox, mass_vapourised = self.mass_evapouration(
                self.tank.volume_tank_oxidizer, mass_liquid_ox_old, self.mass_vapour_ox[iteration], ox_density_vapour, ox_density_liquid)

                
        Heat_of_vapourisation = self.oxidizer.get_enthalpy_vapour(temperature) - self.oxidizer.get_enthalpy_liquid(temperature)
        heat_removed_deltaQ = mass_vapourised * Heat_of_vapourisation
        deltaT = -heat_removed_deltaQ / (mass_liquid_ox_new_it * self.oxidizer.get_cp_liquid(temperature))
        temperature += deltaT  # Change in temperature due to evaporation cooling

        self.dp_dt.append(pressure_old - ox_pressure_liquid)
        self.ox_pressure_liquid.append(ox_pressure_liquid)
        self.ox_liquid_massflowrate.append(ox_liquid_massflowrate)
        self.mass_fuel.append(mass_fuel)
        self.mass_liquid_ox.append(mass_liquid_ox_new_it)
        self.mass_vapour_ox.append(mass_vapour_ox)
        self.temperature_tank.append(temperature)
        self.temperature_combustion.append(self.cea.gas.T)
        self.total_thrust.append(thrust)
        pressure_old = ox_pressure_liquid

        return pressure_old, o_f_old



    def run_vapor_phase_simulation(self, iteration):
        
        self.fuel_pressure.append(0)
        self.fuel_massflowrate.append(0)
        self.mass_fuel.append(0)
        self.mass_liquid_ox.append(0)
        self.O_F_shift.append(0)
        self.pressure_combustion.append(self.config.P_ATMOSPHERE)
        self.ox_liquid_massflowrate.append(0)


        density_vapour = self.mass_vapour_ox[iteration] / (self.tank.volume_tank_oxidizer )

        ox_vapour_massflowrate = self.injector.calculate_mass_flow_rate_SPI( 
            density_vapour, self.ox_pressure_liquid[iteration], self.config.P_EXIT) * self.config.TIME_STEP # nije P_EXIT ni P_ATMOSPHERE zapravo jer se tlak mjenja u kako raketa leti, potrebno integirrati flight dynamics za bolju tocnost

        # mv_dot = N_ORFICE_OX * C_D * A_ox_injector * np.sqrt(2 * density_v * (pressure_v - P_ATMOSPHERE)) * TIME_STEP
        vapour_mass_old = self.mass_vapour_ox[iteration]
        mass_vapour_ox_new = self.mass_vapour_ox[iteration] - ox_vapour_massflowrate

        gamma = self.oxidizer.get_gamma_gas()

        # Update temperature and pressure using polytropic relations
        temperature_new = self.temperature_tank[iteration] * (mass_vapour_ox_new / vapour_mass_old) ** (gamma - 1)
        pressure_new = self.ox_pressure_liquid[iteration] * (temperature_new / self.temperature_tank[iteration]) ** (gamma / (gamma - 1))

        # print("temperatura",temperature_new, "pressure", pressure_new, "mass_vapour_ox_new", mass_vapour_ox_new, "mass_vapour_ox_old", vapour_mass_old)

        self.ox_pressure_liquid.append(pressure_new)
        self.mass_vapour_ox.append(mass_vapour_ox_new)
        self.temperature_tank.append(temperature_new)


    def mass_evapouration(self, tank_volume, mass_liquid_ox_old, mass_vapour_ox, ox_density_vapour, ox_density_liquid):
        # http://www.aspirespace.org.uk/downloads/The%20physics%20of%20nitrous%20oxide.pdf
        mass_liquid_ox = (tank_volume - (mass_liquid_ox_old + mass_vapour_ox) / ox_density_vapour ) / ( 1 / ox_density_liquid - 1 / ox_density_vapour)
        mass_vapourised = mass_liquid_ox_old - mass_liquid_ox
        mass_vapour_ox += mass_vapourised
        return mass_liquid_ox, mass_vapour_ox, mass_vapourised
    
    def calculate_new_combustion_pressure(self, k, combustion_temperature, sumOf_masflowrates, Area_throat, R):
        return sumOf_masflowrates/Area_throat * np.sqrt(combustion_temperature) *np.sqrt(R/k) * ((k+1)/2)**((k + 1) / (2 * (k - 1)))
    
    def update_time(self, iteration):
            time = iteration * self.config.TIME_STEP
            self.time.append(time)

    def calculate_escape_velocity(self, k, Tc, M, pe, pc):
        Rmolar = 8314  # J/(kmol*K)
        return np.sqrt(2 * (Rmolar * k) / (k - 1) * Tc / M * (1 - (pe / pc) ** ((k - 1) / k)))


    def plot_simulation_results(self):
        # Convert lists to numpy arrays for consistency
        fig = 9.5
        plt.figure(figsize=(fig * 2, fig))

        # Plot for tank and chamber pressures
        plt.subplot(2, 2, 1)
        plt.plot(self.time, self.ox_pressure_liquid , color='red', label='Oxidizer Pressure')
        plt.plot(self.time, self.pressure_combustion , color='blue', label='Chamber Pressure')
        plt.plot(self.time, self.fuel_pressure , color='green', label='Fuel Pressure')
        plt.xlabel('Time (s)')
        plt.ylabel('Pressure (bar)')
        plt.legend()
        plt.grid(True)

        # Plot for mass flow rates
        plt.subplot(2, 2, 2)
        ox_liquid_massflowrate_np = np.array(self.ox_liquid_massflowrate)/ self.config.TIME_STEP
        fuel_massflowrate_np = np.array(self.fuel_massflowrate) / self.config.TIME_STEP
        plt.plot(self.time, ox_liquid_massflowrate_np , color='red', label='N2O Mass Flow Rate')
        plt.plot(self.time, fuel_massflowrate_np , color='green', label='Fuel Mass Flow Rate')
        plt.xlabel('Time (s)')
        plt.ylabel('Mass Flow Rate (kg/s)')
        plt.legend()
        plt.grid(True)

        # Plot for mass in tank
        plt.subplot(2, 2, 3)
        plt.plot(self.time, self.mass_liquid_ox, color='red', label='Mass Liquid')
        plt.plot(self.time, self.mass_vapour_ox, color='green', label='Mass Vapour')
        plt.xlabel('Time (s)')
        plt.ylabel('Mass in Tank (kg)')
        plt.legend()
        plt.grid(True)

        # Plot for total thrust
        plt.subplot(2, 2, 4)
        plt.plot(self.time[:len(self.total_thrust)], self.total_thrust, color='blue')
        plt.xlabel('Time (s)')
        plt.ylabel('Total Thrust (N)')
        plt.grid(True)

        plt.tight_layout()
        plt.show()