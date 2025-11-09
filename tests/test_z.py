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
ensure_lib("mpmath")
ensure_lib("matplotlib")

import mpmath as mp
import matplotlib.pyplot as plt

mp.dps = 250 # precision

# Z
def Z(t):
    s = mp.mpc(0.5, t)
    z_val = mp.zeta(s)
    return z_val

t = 11223344.6743864181695953551720
z = Z(t)
abs_z = abs(z)

print(abs_z)