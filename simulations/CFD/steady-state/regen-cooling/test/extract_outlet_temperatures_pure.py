import pyvista as pv
import numpy as np
import os

# --- Settings ---
region = "fluid-domena"
time_step = "1000"
patch_name = "outlet"

# --- Correct VTP path ---
vtp_file = f"VTK/{region}/simpleMlaznica_{time_step}/boundary/{patch_name}.vtp"

if not os.path.exists(vtp_file):
    raise FileNotFoundError(f"Cannot find: {vtp_file}")

print(f"[INFO] Reading: {vtp_file}")
mesh = pv.read(vtp_file)

# --- Extract temperature field ---
if 'T' not in mesh.point_data:
    raise KeyError("Temperature field 'T' not found in patch data.")

T = mesh.point_data['T']

# --- Compute stats ---
print("\n--- Outlet Temperature Stats ---")
print(f"Min  : {np.min(T):.4f}")
print(f"Max  : {np.max(T):.4f}")
print(f"Mean : {np.mean(T):.4f}")
print(f"Count: {len(T)} points")



from CoolProp.CoolProp import PhaseSI, PropsSI

def ipa_phase(T, P):
    try:
        # Attempt to use CoolProp with REFPROP backend if installed
        return PhaseSI('P', P * 1e5, 'T', T, 'REFPROP::IPA')
    except ValueError:
        # Fallback to critical properties if REFPROP isn't available
        T_crit = PropsSI('Tcrit', 'IPA')
        P_crit = PropsSI('pcrit', 'IPA') / 1e5  # Convert Pa to bar
        
        if T > T_crit and P > P_crit:
            return 'supercritical'
        elif T > T_crit:
            return 'supercritical_gas'
        elif P > P_crit:
            return 'supercritical_liquid'
        else:
            return PhaseSI('P', P * 1e5, 'T', T, 'IPA')

# Usage (600 K = 326.85°C, 50 bar)
print(ipa_phase(600, 50))  # Output depends on fluid availability
