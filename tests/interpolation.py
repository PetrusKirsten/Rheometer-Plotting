import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Shear rates from the first dataset (300 to 0.1 1/s) where you want to interpolate shear stress
shear_rate_1 = np.array([
    300.0204,
    278.6028,
    257.1793,
    235.7568,
    214.3328,
    192.9096,
    171.4863,
    150.0632,
    128.6399,
    107.2167,
    85.7933,
    64.3701,
    42.9465,
    21.5235,
    0.1000])

# Shear rates from the second dataset (100 to 0.1 1/s)
shear_rate_2 = np.array([
    100.0105972, 92.87416077, 85.7375946, 78.60144043, 71.46489716,
    64.32920837, 57.19113159, 50.05503464, 42.91851425, 35.78365707,
    28.64353752, 21.51111221, 14.37381363, 7.236902714, 0.069915041])

# Shear stress values for the second dataset
shear_stress_2 = np.array([
300.5299377,
296.9217529,
293.1928711,
288.539093,
284.1900024,
278.8563843,
273.988739,
268.7279663,
262.3976746,
254.9573975,
245.9840698,
233.469635,
217.4181671,
192.7894745,
111.1317749,
])

# Interpolation function (linear interpolation) for the second dataset
interpolation_function = interp1d(shear_rate_2, shear_stress_2, kind='linear', fill_value="extrapolate")

# Interpolate shear stress values for the first dataset's shear rates using the second dataset's data
interpolated_shear_stress = interpolation_function(shear_rate_1)

# Output interpolated shear stress results
print("Interpolated Shear Stress Results (from 2nd dataset) for the First Dataset's Shear Rates:")
for sr, ss in zip(shear_rate_1, interpolated_shear_stress):
    print(f"{ss}")
    # print(f"Shear rate: {sr:.4f} 1/s, Estimated shear stress: {ss:.4f} Pa")

# Output interpolated shear stress results
print("Interpolated Shear Stress Results (from 2nd dataset) for the First Dataset's Shear Rates:")
for sr, ss in zip(shear_rate_1, interpolated_shear_stress):
    print(f"{sr}")
    # print(f"Shear rate: {sr:.4f} 1/s, Estimated shear stress: {ss:.4f} Pa")

# Optional: Plotting the original second dataset and the interpolated results
plt.figure(figsize=(10, 6))
plt.plot(shear_rate_2, shear_stress_2, 'o', label="Second Dataset (Original)")
plt.plot(shear_rate_1, interpolated_shear_stress, 'o', label="Interpolated for First Dataset Shear Rates")
plt.xlabel("Shear Rate (1/s)")
plt.ylabel("Shear Stress (Pa)")
plt.title("Interpolated Shear Stress for First Dataset's Shear Rates (from Second Dataset)")
plt.legend()
plt.xscale('linear')
plt.grid(True)
plt.show()
