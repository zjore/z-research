#!/usr/bin/env python3
import subprocess
import sys

def ensure_lib(package):
    try:
        __import__(package)
    except ImportError:
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", package, "--quiet"]
        )

# Ensure required libraries
ensure_lib("numpy")
ensure_lib("matplotlib")
ensure_lib("pandas")

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

# ----------------- Load zeros (imag parts) -----------------
df = pd.read_csv("datasets/5000_zeros_for_testing.csv")

# Expect a column named 'real_y' with gamma values > 0
y_values = pd.to_numeric(df["real_y"], errors="coerce").to_numpy(dtype=np.float64)
y_values = y_values[np.isfinite(y_values)]
y_values = y_values[y_values > 0.0]
y_values.sort()

# Optional: cap how many zeros to include
# N = 5000
# y_values = y_values[:N]

print(f"Using {len(y_values)} zeros")

# ----------------- Grid and correction -----------------
x = np.linspace(1.0, 50.0, 1000, dtype=np.float64)
logx = np.log(x)

# Initialize as complex128
correction = np.zeros_like(x, dtype=np.complex128)

# Loop over zeros (kept for simplicity)
for y in y_values:
    rho = 0.5 + 1j * y
    # numerically stable complex power: exp(rho * log x)
    term = np.exp(rho * logx) / rho
    correction += term

# Real part with conjugate doubling
correction_real = -2.0 * correction.real

# ----------------- Plot -----------------
plt.figure(figsize=(12, 6))
plt.plot(x, x, label="Main Term (x)", linestyle="dashed", color="blue")
plt.plot(x, x + correction_real, label=f"Main Term + Correction (N={len(y_values)})", color="red")

plt.xlabel("x")
plt.ylabel("Approximation of ψ(x)")
plt.title("Approximation of Prime Counting with Fine Grid")

plt.minorticks_on()
plt.grid(which="major", linestyle="-", linewidth=0.5, color="grey", alpha=0.7)
plt.grid(which="minor", linestyle=":", linewidth=0.5, color="lightgrey", alpha=0.5)

plt.xticks(np.arange(0, np.max(x) + 1, 5))
plt.yticks(np.arange(0, np.max(x) + 1, 5))

plt.legend()
plt.tight_layout()

plt.show()
