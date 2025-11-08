#!/usr/bin/env python3
import subprocess, sys
def ensure_lib(package):
    try: __import__(package)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package, "--quiet"])

# deps
for p in ("numpy","matplotlib","pandas"): ensure_lib(p)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ----------------- Load zeros (imag parts γ_n) -----------------
df = pd.read_csv("datasets/5000_zeros_for_testing.csv")
gammas = pd.to_numeric(df["real_y"], errors="coerce").to_numpy(dtype=np.float64)
gammas = gammas[np.isfinite(gammas)]
gammas = gammas[gammas > 0.0]
gammas.sort()
print(f"Using {len(gammas)} zeros")

# ----------------- Grid -----------------
x = np.linspace(1.0, 50.0, 1000, dtype=np.float64)
logx = np.log(x)
sqrtx = np.sqrt(x)

# ----------------- Real-form correction (no complex math) -----------------
# Each γ contributes:  - x^{1/2} * [cos(γ log x) + 2γ sin(γ log x)] / (γ^2 + 1/4)
def correction_real_direct(gammas, x, logx, sqrtx, chunk_size=2000):
    corr = np.zeros_like(x, dtype=np.float64)
    n = len(gammas)
    for i in range(0, n, chunk_size):
        g = gammas[i:i+chunk_size]                      # shape (G,)
        D = (g*g + 0.25)[:, None]                       # (G,1)
        theta = g[:, None] * logx[None, :]              # (G,X)
        # real pair contribution (with explicit minus sign for the explicit formula)
        block = - sqrtx[None, :] * (np.cos(theta) + 2.0*g[:, None]*np.sin(theta)) / D
        corr += block.sum(axis=0)
    return corr

correction_real = correction_real_direct(gammas, x, logx, sqrtx, chunk_size=2000)

# ----------------- Plot -----------------
plt.figure(figsize=(12, 6))
plt.plot(x, x, label="Main Term (x)", linestyle="dashed")
plt.plot(x, x + correction_real, label=f"Main Term + Correction (N={len(gammas)})")

plt.xlabel("x")
plt.ylabel(r"Approximation of $\psi(x)$")
plt.title("Prime-Counting Approximation via Conjugate-Paired Zero Sum (Real Form)")

plt.minorticks_on()
plt.grid(which="major", linestyle="-", linewidth=0.5, alpha=0.7)
plt.grid(which="minor", linestyle=":", linewidth=0.5, alpha=0.5)

import math
maxv = math.ceil(x.max())
plt.xticks(np.arange(0, maxv+1, 5))
plt.yticks(np.arange(0, maxv+1, 5))

plt.legend()
plt.tight_layout()
plt.show()
