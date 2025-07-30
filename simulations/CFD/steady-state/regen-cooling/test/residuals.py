import re
import matplotlib.pyplot as plt

# File path
log_file = 'log.chtMultiRegionSimpleFoam'

# Patterns
time_pattern = re.compile(r'^Time = ([\d\.]+)')
fluid_start_pattern = re.compile(r'^Solving for fluid region')
residual_pattern = re.compile(r'Solving for (\w+), Initial residual = ([\deE\+\.-]+)')

# Data structures
residuals = {'Ux': [], 'Uy': [], 'Uz': [], 'h': [], 'p_rgh': []}
time_values = []
current_time = None
parsing_fluid = False
temp_residuals = {}

# Read and parse
with open(log_file, 'r') as f:
    for line in f:
        # Time step
        time_match = time_pattern.search(line)
        if time_match:
            if temp_residuals:  # Store previous time step data
                time_values.append(current_time)
                for key in residuals:
                    residuals[key].append(temp_residuals.get(key, None))
                temp_residuals.clear()
            current_time = float(time_match.group(1))
            parsing_fluid = False  # Reset region context

        # Detect start of fluid region block
        if fluid_start_pattern.search(line):
            parsing_fluid = True

        # Extract residuals in fluid block only
        if parsing_fluid:
            res_match = residual_pattern.search(line)
            if res_match:
                var, res = res_match.groups()
                if var in residuals:
                    temp_residuals[var] = float(res)

# Append last timestep
if current_time and temp_residuals:
    time_values.append(current_time)
    for key in residuals:
        residuals[key].append(temp_residuals.get(key, None))

# Plotting
plt.figure(figsize=(10, 6))
for var, values in residuals.items():
    # Remove None values for safe plotting
    clean_time = [t for t, v in zip(time_values, values) if v is not None]
    clean_vals = [v for v in values if v is not None]
    plt.semilogy(clean_time, clean_vals, label=var)

plt.xlabel('Time')
plt.ylabel('Initial Residual (log scale)')
plt.title('Fluid Region Initial Residuals from log.chtMultiRegionSimpleFoam')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
