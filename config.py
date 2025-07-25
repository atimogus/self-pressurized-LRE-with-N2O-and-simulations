class Config:
    # CEA Configuration
    YAML_FILE = "gri30.yaml"
    FUEL = "Acetone"        # closest fuel to IPA in coolprop library
    OXIDIZER = "N2O"
    FUEL_FOR_CEA = 'C3H8'   # closest fuel to IPA in cantera library 

    # Simulation Parameters
    TIME_STEP = 0.03
    INITIAL_TEMP_DEGREES = 23                        # degrees Celsius
    INITIAL_TEMP_K = INITIAL_TEMP_DEGREES + 273.15  # Kelvin
    gravity_constant = 9.81 # m/s^2, gravitational constant

    # Tank Properties
    OXIDIZER_MASS_KG = 22   
    FUEL_DURATION = 1.1    # how much shorter will fuel last than oxidizer (percentage of oxidizer duration) 0.5 is half of oxidizer duration
    OXIDIZER_TANK_DIAMETER = 160
    OXIDIZER_TANK_THICKNESS = 5 #mm
    FUEL_TANK_DIAMETER = 70 #milimeters
    FUEL_TANK_THICKNESS = 2 #mm
    TANK_ULAGE = 0.05 # percentage of tank volume that is not filled with liquid (filled with vapour)
    PRESSURE_FRICTION_LOSS = 2e5  # pascals


    # Injector Properties
    DISCHARGE_COEFFICIENT = 0.66  # discharge coefficient
    discharge_coefficient_fuel = 0.66
    OXIDIZER_ORIFICE_DIAMETER = 2  # mm
    NUMBER_OF_ORIFICES_OXIDIZER = 9
    NUMBER_OF_ORIFICES_FUEL = 12
    # PINTLE_ANGLE = 60

    # Chamber Properties
    OXIDIZER_TO_FUEL_RATIO = 3.2
    COMBUSTION_CHAMBER_PRESSURE = 28e5  # pascals (pressure of combustion in chamber)
    CHAMBERSAFE_THICKNESS = 5  # mm

    # Exhaust Properties
    P_ATMOSPHERE = 101325  # pascals
    P_EXIT = 95555          #pascals

    # Flight Parameters
    DRY_WEIGHT = 25  # kilograms (without payload)
    LAUNCH_RAIL_ANGLE = 84  # degrees
    DRAG_COEFFICIENT = 0.7  # for mach M, Cd is = https://www.nakka-rocketry.net/RD_body.html
    LOSS = 7  # percent
    L_D_RATIO = 8.5  # tank length / body diameter ratio


    CONVERGENT_NOZZLE_ANGLE = 45  # degrees - dio koj se smanjuje
    DIVERGENT_NOZZLE_ANGLE = 16.6  # degrees - dio koj se povecava
    combustion_chamber_height = 150 # mm

    TANK_MATERIAL_STRENGTH = 230e6
    SAFETY_FACTOR_TANK = 2
    SAFETY_FACTOR_CHAMBER = 3
    # Fuel Gas Orifice Sizing
  
    GAS_VELOCITY = 1  # m/s -> this should not be higher than Mach=1

    NUNMER_OF_REGEN_CHANNELS = 12
    ALU_THICKNESS = 2e-3  # m
    D_HIDRAULIC_HOSE = 25  # mm

    execute_simulation = False  # if True, simulation that is defined below will be runned, if False no neither coldflow or hot fire test will run
    Cold_Flow_test = False # if True, run the cold flow test if false, run the hot fire test
    flight_simulation = False # if True, run the flight simulation and Cold_Flow_test must be false so Hot fire test can run

    generate_Freecad_geometry = False

    nozzle_contour_optimization = True
    num_of_points_nozzle = 100