"""
ME424_Python_Design_study.py 
Output:  printed table and design_study__results.csv
Needs:   numpy and the two design-day CSV files in the same folder.
"""

import numpy as np

# ----------------------------------------------------------------------------
# 1. FIXED INPUTS
# ----------------------------------------------------------------------------
A        = 8.0      # collector area, m2
V_TANK   = 0.150    # tank volume, m3
RHO      = 1000.0   # water density, kg/m3
CP       = 4180.0   # water specific heat, J/kgK
MU       = 5.77e-4  # water viscosity, Pa s (for pump pressure drop)
K_W      = 0.64     # water conductivity, W/mK
T_START  = 10.0     # mains water temperature, C
T_TARGET = 45.0     # target temperature, C
DT       = 60.0     # time step, s (1 minute)

ALPHA, TAU    = 0.90, 1.00     # absorptance, transmittance (brief)
H_SLOPE       = 1.757          # tube length = module length up the slope of roof, m
PACKING       = 0.996          # fraction of the module covered by cells
ETA_REF, BETA = 0.243, 0.0026  # PV efficiency at 25 C, and its drop per K (datasheet)

# Heat loss coefficients from the top of the collector (all W/m2K)
H_RAD   = 5.5            # radiation to the sky: 4*eps*sigma*T^3 at about 300 K
H_WIND  = 15.0           # wind convection on the outer surface
H_GAP   = 3.0            # convection across the air gap under the glass (Hollands, 25 mm gap)
R_GLASS = 0.002 / 1.0    # conduction resistance of 2 mm glass (thickness / k), m2K/W

# Baseline design. Each test copies this and changes one value.
BASE = dict(glazed=True, W=0.10, D=0.010, mdot=0.05)   # glass?, tube spacing m, tube diameter m, flow kg/s


# ----------------------------------------------------------------------------
# 2. WEATHER
# ----------------------------------------------------------------------------
# Read a design day (hour, sun W/m2, air temp C) and spread it onto a 1-minute grid.
def load_day(filename):
    data = np.loadtxt(filename, delimiter=",", comments="#")
    hours = np.arange(5.0, 20.0, DT / 3600) 
    return hours, np.interp(hours, data[:, 0], data[:, 1]), np.interp(hours, data[:, 0], data[:, 2])


# ----------------------------------------------------------------------------
# 3. COLLECTOR MODEL (one moment in time, all temperatures in C)
# ----------------------------------------------------------------------------
def sky_temp(Ta):
    return 0.0552 * (Ta + 273.15) ** 1.5 - 273.15     # Swinbank sky temperature, C


def top_loss(d):
    # Top loss coefficient: gap, glass and wind in series if glazed, otherwise just wind.
    return 1 / (1 / H_GAP + R_GLASS + 1 / H_WIND) if d["glazed"] else H_WIND


def collector(G, Ta, T_in, mdot, d):
    """Returns heat to water (W), electricity (W) and plate temperature (C)."""
    T_sky = sky_temp(Ta)
    U_top = top_loss(d)

    # Plate-to-water: laminar h = 3.66 k / D, over tube area pi*D per spacing W.
    # D cancels, so only the tube spacing matters.
    U_w = 3.66 * K_W * np.pi / d["W"]               # W/m2K of plate
    if mdot > 0:
        NTU = U_w * A / (mdot * CP)                 # size of collector vs flow heat capacity
        k = mdot * CP / A * (1 - np.exp(-NTU))      # effectiveness method, W/m2K
    else:
        k = 0.0                                     # pump off: no heat to water

    S = TAU * ALPHA * G                 # sun absorbed, W/m2
    c = PACKING * ETA_REF * TAU * G     # electricity, W/m2, if the plate were at 25 C

    # Plate energy balance (energy in = energy out), every term is linear in Tp:
    #   S - c(1 - BETA(Tp - 25)) - H_RAD(Tp - T_sky) - U_top(Tp - Ta) - k(Tp - T_in) = 0
    # Rearranged to give Tp directly:
    Tp = ((S - c * (1 + 25 * BETA) + H_RAD * T_sky + U_top * Ta + k * T_in)
          / (H_RAD + U_top + k - c * BETA))

    q = k * (Tp - T_in) * A                  # heat to water, W
    elec = c * (1 - BETA * (Tp - 25)) * A    # electricity, W (falls as the plate heats up)
    return q, elec, Tp


# Pump electrical power from the pressure drop (Darcy-Weisbach), W.
def pump_power(d, mdot, L_pipe=20.0, D_pipe=0.0136, K_fit=10.0, eta_pump=0.25):
    def dp(m, D, L, K):
        v = m / (RHO * np.pi * D**2 / 4)                  # velocity, m/s
        Re = 4 * m / (np.pi * D * MU)                     # Reynolds number
        f = 64 / Re if Re < 2300 else 0.316 * Re**-0.25   # laminar or Blasius friction factor
        return (f * L / D + K) * RHO * v**2 / 2           # pipe friction + fittings, Pa
    n_tubes = A / (d["W"] * H_SLOPE)
    total = dp(mdot / n_tubes, d["D"], H_SLOPE, 1.5) + dp(mdot, D_pipe, L_pipe, K_fit)   # one tube + pipework
    return total * (mdot / RHO) / eta_pump               # pressure x flow / efficiency


# ----------------------------------------------------------------------------
# 4. STRATEGY A: recirculate the whole tank, pump on only when it gains heat
# ----------------------------------------------------------------------------



# ----------------------------------------------------------------------------
# 5. STRATEGY B (estimate): fill the tank once with 45 C water
# ----------------------------------------------------------------------------



# ----------------------------------------------------------------------------
# 6. RUN THE STUDY
# ----------------------------------------------------------------------------