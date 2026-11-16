"""
ME424_Python_Design_study_Final_WORKING_model.py
ID:      23405902
Name:    Oisin Naughton
Output:  printed table and design_study__results.csv
Needs:   the two design-day CSV files in the same folder.
"""

import os, sys

# Work in this script's own folder, so the CSV files are found even when the script is
# double-clicked (Windows may otherwise start Python in C:\Windows\System32).
os.chdir(os.path.dirname(os.path.abspath(__file__)))


# If anything goes wrong, show the error and wait, instead of the window closing instantly.
def pause_on_error(*error):
    sys.__excepthook__(*error)
    input("\nError above. Press Enter to close...")
sys.excepthook = pause_on_error


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
def run_recirc(d, day):
    """Returns tank end temp (C), time 45 C reached, hours from pump start, electricity kWh."""
    hours, G, Ta = day
    T, elec, t_on, t45 = T_START, 0.0, None, None
    for h, g, ta in zip(hours, G, Ta):                 # step through the day minute by minute
        q, e, Tp = collector(g, ta, T, d["mdot"], d)
        if q <= 0:                                     # collector colder than tank: pump off
            q, e, Tp = collector(g, ta, T, 0.0, d)     # (q is 0 with no flow)
        elif t_on is None:
            t_on = h                                   # first time the pump runs
        T += q * DT / (V_TANK * RHO * CP)              # warm the mixed tank: dT = Q / (m cp)
        elec += e * DT
        if t45 is None and T >= T_TARGET:
            t45 = h                                    # first time the tank reaches 45 C
    return T, t45, (None if t45 is None else t45 - t_on), elec / 3.6e6   # J -> kWh


# ----------------------------------------------------------------------------
# 5. STRATEGY B (estimate): fill the tank once with 45 C water
# ----------------------------------------------------------------------------
def run_oncethrough(d, day):
    """The flow is throttled so water leaves at 45 C. With such a slow flow the plate sits at
    about 45 C, so useful heat = sun absorbed - electricity - losses, all with Tp = 45 C.
    Water made = Q / (cp * 35 K). Returns litres made and time the tank is full."""
    hours, G, Ta = day
    Tp = T_TARGET                                      # plate at about 45 C
    litres, t_full = 0.0, None
    for h, g, ta in zip(hours, G, Ta):
        q = (TAU * ALPHA * g                                   # sun absorbed
             - PACKING * ETA_REF * TAU * g * (1 - BETA * (Tp - 25))   # electricity
             - H_RAD * (Tp - sky_temp(ta))                     # radiation to sky
             - top_loss(d) * (Tp - ta)) * A                    # convection to air, W total
        if q > 0:                                      # plate can get hotter than 45 C
            litres += q * DT / (CP * (T_TARGET - T_START))   # kg of water = litres
        if t_full is None and litres >= V_TANK * 1000:
            t_full = h
    return litres, t_full


# Print a decimal hour (13.5) as a clock time ("13:30").
def clock(h):
    return "not reached" if h is None else f"{int(h):02d}:{int(round((h % 1) * 60)) % 60:02d}"


# ----------------------------------------------------------------------------
# 6. RUN THE STUDY
# ----------------------------------------------------------------------------
DAYS = {"Average": load_day("design_day_average_south35.csv"),
        "Clear":   load_day("design_day_clear_south35.csv")}
rows = []   # results for the CSV file


# Run one design on both days, print one table row and save it.
def report(test, value, d):
    line = f"{test:<12}{value:<10}"
    for name, day in DAYS.items():
        T_end, t45, hours, kwh_e = run_recirc(d, day)
        line += f"{T_end:>8.1f}{clock(t45):>13}{kwh_e:>7.1f}"
        rows.append([test, value, name, round(T_end, 1), clock(t45),
                     "-" if hours is None else round(hours, 1), round(kwh_e, 2)])
    print(line + f"{pump_power(d, d['mdot']):>9.1f}")


print(f"{'':22}{'------ Average day ------':>28}{'------- Clear day -------':>28}")
print(f"{'Test':<12}{'Value':<10}" + f"{'T end':>8}{'45 C at':>13}{'kWh e':>7}" * 2 + f"{'Pump W':>9}")

# Test 1: glass cover
for g in (True, False):
    report("Glazing", "glazed" if g else "unglazed", {**BASE, "glazed": g})

# Test 2: flow rate
for m in (0.005, 0.01, 0.02, 0.05, 0.1, 0.2):
    report("Flow kg/s", m, {**BASE, "mdot": m})

# Test 3: tube spacing
for W in (0.05, 0.10, 0.15):
    report("Spacing mm", int(W * 1000), {**BASE, "W": W})

# Test 4: tube diameter (laminar flow, so heat transfer is the same; only pump power changes)
print("\nTube diameter (only pump power changes):")
for D in (0.006, 0.008, 0.010, 0.015):
    d = {**BASE, "D": D}
    n_tubes = A / (d["W"] * H_SLOPE)
    Re = 4 * (d["mdot"] / n_tubes) / (np.pi * D * MU)      # below 2300 = laminar
    print(f"  D = {D*1000:4.0f} mm   Re = {Re:5.0f}   pump power = {pump_power(d, d['mdot']):5.2f} W")

# Test 5: circulation strategy
print("\nCirculation strategy (baseline design):")
for name, day in DAYS.items():
    T_end, t45, _, _ = run_recirc(BASE, day)
    litres, t_full = run_oncethrough(BASE, day)
    print(f"  {name:8} A recirculate: tank {T_end:5.1f} C, 45 C at {clock(t45)}"
          f"  |  B once-through (estimate): {litres:5.0f} L at 45 C, tank full at {clock(t_full)}")

# Save the table to a CSV file.
with open("design_study_results.csv", "w") as f:
    f.write("test,value,day,tank_end_C,time_45C_UTC,hours_to_45C,electricity_kWh\n")
    for r in rows:
        f.write(",".join(str(x) for x in r) + "\n")
print("\nSaved design_study_results.csv")

# Keep the window open when the script is double-clicked.
input("\nPress Enter to close...")
