from rocketpy import Fluid, MassBasedTank, Environment, Flight, Rocket, LiquidMotor

class rocketpy_sim:
    #https://docs.rocketpy.org/en/latest/reference/index.html
    #https://colab.research.google.com/drive/1Pzf_dogoRFYmhslLHxZb6HzG9On6ghZh
     def __init__(self, oxidizer, fuel, config, cea, tank, injector, nozzle):

    #define fluids
    temp, density_v, pressure_v, density_l, pressure_l = calculate_initial_properties(INITIAL_TEMP_DEGREES, OXIDIZER, GAS)
    density_f = PropsSI('D', 'P', pressure_l, 'T', temp, FUEL)

    oxidizer_liq = Fluid(
        name="Liquid N2O",
        density=density_l # kg/m3
        )
    oxidizer_gas = Fluid(
        name="Gas N2O",
        density=density_v # kg/m3, program estimates tank overfill if number is lower than 5, i put 10 for safety reasons
    )
    fuel_liq = Fluid(
        name="Liquid ethanol",
        density=density_f # kg/m3
        )
    #define tank
    ox_tanks_geometry = CylindricalTank(
        radius=body_diameter/2, # meters (or  3.94 in)
        height=tank_height, # meters (or 47.24 in)
        spherical_caps=False
    )
    
    fuel_tanks_geometry = CylindricalTank(
        radius=fuel_tank_diameter/2, # meters (or  3.94 in)
        height=tank_height, # meters (or 47.24 in)
        spherical_caps=False
    )
    
    V_liquid, V_gas, V_tank, mass_v, mass_l = calculate_tank_properties(density_l, density_v)
    
    tank_ox_mass_array = np.array(data["mass_liquid"])
    time_data = np.arange(0, len(tank_ox_mass_array) * TIME_STEP, TIME_STEP)
    liquid_massflow_2d = [[time, flow] for time, flow in zip(time_data, data["mass_liquid"])]
    
    tank_oxgas_mass_array = np.array(data["mass_vapour"])
    time_data = np.arange(0, len(tank_oxgas_mass_array) * TIME_STEP, TIME_STEP)
    gas_ox_massflow_2d = [[time, flow] for time, flow in zip(time_data, data["mass_vapour"])] 
    
    oxidizer_tank = MassBasedTank(name = "oxidizer tank", 
        geometry = ox_tanks_geometry, 
        flux_time = 10, 
        liquid = oxidizer_liq, 
        gas = oxidizer_gas, 
        liquid_mass = liquid_massflow_2d, 
        gas_mass = gas_ox_massflow_2d
    )
    
    # tank_massflow_array = np.array(data["liquid_massflow"])
    # time_data = np.arange(0, len(tank_massflow_array) * TIME_STEP, TIME_STEP)
    # liquid_massflow_2d = [[time, flow] for time, flow in zip(time_data, data["liquid_massflow"])]
    
    # oxidizer_tank = rp.MassFlowRateBasedTank(
    #     name="oxidizer tank",
    #     geometry=ox_tanks_geometry,
    #     flux_time=8, # seconds
    #     initial_liquid_mass=OXIDIZER_MASS_L_KG, # kg
    #     initial_gas_mass=mass_v, # kg
    #     liquid_mass_flow_rate_in=0,
    #     liquid_mass_flow_rate_out= liquid_massflow_2d, # <- defined in the previous cell
    #     gas_mass_flow_rate_in=0,
    #     gas_mass_flow_rate_out=0,
    #     liquid=oxidizer_liq,
    #     gas=oxidizer_gas,
    # )
    # oxidizer_tank.draw()
    # oxidizer_tank.center_of_mass()

    fuel_tank_massflow_array = np.array(data["fuel_massflow"])
    time_data = np.arange(0, len(fuel_tank_massflow_array) * TIME_STEP, TIME_STEP)
    fuel_massflow_2d = [[time, flow] for time, flow in zip(time_data, data["fuel_massflow"])]
    mass_f, V_fuel = calculate_fuel_properties( pressure_l, temp)
    # print(mass_f, mass_v)
    
    fuel_tank = MassBasedTank(name = "fuel tank", 
        geometry = fuel_tanks_geometry, 
        flux_time = 10, 
        liquid = fuel_liq, 
        gas = oxidizer_gas, 
        liquid_mass = fuel_massflow_2d, 
        gas_mass = 0
    )
    # fuel_tank = rp.MassFlowRateBasedTank(
    #     name="fuel tank",
    #     geometry=fuel_tanks_geometry,
    #     flux_time=8, # seconds
    #     initial_liquid_mass=mass_f, # kg
    #     initial_gas_mass=mass_v, # kg
    #     liquid_mass_flow_rate_in=0,
    #     liquid_mass_flow_rate_out=fuel_massflow_2d,
    #     gas_mass_flow_rate_in=0,
    #     gas_mass_flow_rate_out=0,
    #     liquid=fuel_liq,
    #     gas=fuel_gas,
    # )
    # fuel_tank.draw()
    # fuel_tank.center_of_mass()
    
    thrust_array = np.array(data["total_thrust"])
    time_data = np.arange(0, len(thrust_array) * TIME_STEP, TIME_STEP)
    thrust_2d = [[time, flow] for time, flow in zip(time_data, data["total_thrust"])]

    #definiranje motora
    example_liquid_motor = LiquidMotor(
        thrust_source=thrust_2d,
        dry_mass=DRY_WEIGHT*0.085, # kg
        dry_inertia=(0.002, 0.002, 0.002), # kg/m2 , ne mjenja puno stvari
        nozzle_radius=De/(2 * 1000),
        center_of_dry_mass_position=1,
        nozzle_position=0,
    )
    
    # example_liquid_motor = SolidMotor(
    #     thrust_source = thrust_2d,
    #     dry_mass=1,
    #     dry_inertia=(0.0000000000001, 0.0000000000001, 0.0000000000001), 
    #     nozzle_radius=De/(2 * 1000),
    #     grain_number = 1, 
    #     grain_density = 1000, 
    #     grain_outer_radius = De/(2 * 1000), 
    #     grain_initial_inner_radius = De/(4 * 1000), 
    #     grain_initial_height = 0.3, 
    #     grain_separation=0.006, 
    #     grains_center_of_mass_position=-0.683, 
    #     center_of_dry_mass_position=-0.683
    # )

    # example_liquid_motor = GenericMotor(
    #     thrust_source = thrust_2d, 
    #     burn_time = 4,
    #     chamber_radius = Dc_mm /(2* 1000), 
    #     chamber_height = l_chamber, 
    #     chamber_position = 0.1, 
    #     propellant_initial_mass = OXIDIZER_MASS_L_KG + mass_f, 
    #     nozzle_radius = De/(2 * 1000) 
    # )


    # Add the tanks to the motor!
    
    example_liquid_motor.add_tank(
        tank=oxidizer_tank,
        position=0.5 + tank_height/2, # meters
    )
    
    example_liquid_motor.add_tank(
        tank=fuel_tank,
        position=0.5 + tank_height/2, # meters
    )
    
    # example_liquid_motor.info()
    example_liquid_motor.draw()
    # example_liquid_motor.center_of_mass()
    
    #definiranje uvjeta lansiranja
    env = Environment(
        # kordinate za Zagreb
        latitude=45.81,
        longitude=15.98,
    )
    
    env.set_atmospheric_model(
        type="custom_atmosphere",
        wind_v=4, # m/s
        wind_u=3, # m/s
    )
    
    # env.info()
    
    #crtanje ostatka tijela rakete
    example_rocket = Rocket(
        radius=body_diameter/2,
        mass=DRY_WEIGHT,
        inertia=(25, 25, 1),
        power_off_drag=1,
        power_on_drag=1,
        center_of_mass_without_motor=0.5 + tank_height,
        coordinate_system_orientation="nose_to_tail",
    )
    
    example_rocket.add_motor(example_liquid_motor, position= 2 * tank_height)
    example_rocket.add_nose(length=0.7, kind="vonKarman", position=0)  
    example_rocket.add_trapezoidal_fins(
        n=4,
        root_chord=0.355,
        tip_chord=0.0803,
        span=0.156,
        position=2 * tank_height,
    )
    example_rocket.add_parachute(
        name="Drogue",
        cd_s=5,
        trigger="apogee",
    )
    
    example_rocket.info()
    example_rocket.draw()
    
    test_flight = Flight(
        environment=env,
        rocket=example_rocket,
        rail_length=10,
        heading=0,
        inclination=84
    )
    
    test_flight.prints.burn_out_conditions()
    test_flight.prints.apogee_conditions()
    test_flight.prints.impact_conditions()
    test_flight.prints.maximum_values()
    test_flight.plots.trajectory_3d()
    