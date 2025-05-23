from CoolProp.CoolProp import PropsSI

def get_ethanol_critical_point():
    """Returns critical temperature (K) and pressure (bar) for ethanol"""
    try:
        T_crit = PropsSI('Tcrit', 'Ethanol')  # Critical temperature in K
        P_crit = PropsSI('Pcrit', 'Ethanol') / 1e5  # Convert Pa to bar
        return T_crit, P_crit
    except ImportError:
        print("Error: CoolProp not installed. Run 'pip install CoolProp' first.")
        return None, None

# Get and display values
T_max, P_max = get_ethanol_critical_point()
if T_max and P_max:
    print(f"Ethanol critical point:")
    print(f"- Maximum temperature for liquid/vapor: {T_max:.2f} K ({T_max-273.15:.2f}°C)")
    print(f"- Corresponding pressure: {P_max:.2f} bar")
