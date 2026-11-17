"""
ME424_pvt_comparison.py
ID: 23405902
Name: Oisín Naughton
Output:  printed table, comparison_results.csv and comparison_plot.png
Needs:   numpy, matplotlib, the two design-day CSVs and design_study_results.csv

Benchmark comparison (report Section 6): the chosen PV/T design against
    - plain PV      (same 4 Aiko modules, cell temperature from NOCT)
    - plain thermal (Viessmann Vitosol 200-FM, certified efficiency curve)
      on the same 8 m2, 150 L tank and September design days.

The PV/T results are NOT recalculated here: they are read from the baseline
row of design_study_results.csv, so run ME424_Python_Design_study_Final_WORKING_model.py first.
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

import csv
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# 1. INPUTS
# ----------------------------------------------------------------------------
A      = 8.0        # collector area, m2
V_TANK = 0.150      # tank volume, m3
RHO      = 1000.0   # water density, kg/m3
CP       = 4180.0   # water specific heat, J/kgK
T_START  = 10.0     # mains water temperature, C
T_TARGET = 45.0     # target temperature, C
MDOT   = 0.05       # flow rate for the thermal collector, kg/s (same as the PV/T design)
DT     = 60.0       # time step, s

# Plain PV: Aiko Neostar 2P A485 datasheet
PACKING = 0.996     # fraction of the module covered by cells
ETA_REF, BETA = 0.243, 0.0026  # PV efficiency at 25 C, and its drop per K (datasheet)
NOCT    = 25 + (1 - 365 / (485 * 0.8)) / BETA   # from 365 W at NOCT vs 485 W at STC, about 48 C

# Plain thermal: Viessmann Vitosol 200-FM datasheet, aperture area
ETA0, A1, A2 = 0.82, 4.75, 0.024
APERTURE_FRACTION = 0.93   # aperture / gross area (about 2.33 / 2.51 m2)


# ----------------------------------------------------------------------------
# 2. WEATHER (times in UTC)
# ----------------------------------------------------------------------------
# Read a design day (hour, sun W/m2, air temp C) and spread it onto a 1-minute grid.
def load_day(filename):
    data = np.loadtxt(filename, delimiter=",", comments="#")
    hours = np.arange(5.0, 20.0, DT / 3600)
    return hours, np.interp(hours, data[:, 0], data[:, 1]), np.interp(hours, data[:, 0], data[:, 2])


# ----------------------------------------------------------------------------
# 3. PLAIN PV: electricity over the day
# ----------------------------------------------------------------------------
def run_pv(hours, G, Ta):
    T_cell = Ta + (NOCT - 20) / 800 * G                       # cell temperature, C
    power = PACKING * ETA_REF * (1 - BETA * (T_cell - 25)) * G * A   # W at each minute
    return power.sum() * DT / 3.6e6                           # kWh


# ----------------------------------------------------------------------------
# 4. PLAIN THERMAL: heats the tank by recirculation, like the PV/T design
# ----------------------------------------------------------------------------
def run_thermal(hours, G, Ta):
    # Simulates the plain solar-thermal collector heating the tank over one day
    T, heat, t45 = T_START, 0.0, None              # tank starts at 10 C, no heat yet, 45 C not reached
    for h, g, ta in zip(hours, G, Ta):             # step through the day one minute at a time
        q = 0.0                                    # heat this minute (stays 0 if no sun)
        if g > 0:                                  # only calculate when the sun is up
            Tm = T                                 # first guess of mean water temp = tank temp
            for _ in range(5):                     # repeat to refine the mean water temp
                eta = ETA0 - A1 * (Tm - ta) / g - A2 * (Tm - ta)**2 / g   # collector efficiency curve
                q = max(eta, 0) * g * A * APERTURE_FRACTION   # heat into the water (W), aperture area only
                Tm = T + q / (2 * MDOT * CP)       # mean water temp = inlet + half the temp rise
        T += q * DT / (V_TANK * RHO * CP)          # tank temperature rise this minute
        heat += q * DT                             # add this minute's heat to the day total (J)
        if t45 is None and T >= T_TARGET:          # first time the tank reaches 45 C
            t45 = h                                # record that time
    return heat / 3.6e6, T, t45                    # day's heat (kWh), final tank temp, time 45 C reached

def clock(h):
    # Converts a decimal hour (e.g. 12.05) into clock time (e.g. "12:03")
    return "not reached" if h is None else f"{int(h):02d}:{int(round((h % 1) * 60)) % 60:02d}"


# ----------------------------------------------------------------------------
# 5. PV/T RESULTS: baseline row from the design study
# ----------------------------------------------------------------------------
if not os.path.exists("design_study_results.csv"):           # check the design study results file exists
    raise FileNotFoundError("design_study_results.csv not found - run ME424_Python_Design_study.py first.")
pvt = {}                                                     # will hold the PV/T results for each day
with open("design_study_results.csv") as f:                  # open the design study results
    for row in csv.DictReader(f):                            # read it row by row, using the column names
        if row["test"] == "Glazing" and row["value"] == "glazed":   # the baseline design row
            pvt[row["day"]] = row                            # store it under "Average" or "Clear"




# ----------------------------------------------------------------------------
# 6. COMPARE ON BOTH DAYS
# ----------------------------------------------------------------------------
DAYS = {"Average": "design_day_average_south35.csv", "Clear": "design_day_clear_south35.csv"}   # the two design days
rows = []                                                    # results table, one row per system per day
for day, file in DAYS.items():                               # loop over the average and clear days
    hours, G, Ta = load_day(file)                            # read that day's sunshine and air temperature
    p = pvt[day]                                             # PV/T results for that day
    rows.append([day, "PV/T", round(V_TANK * RHO * CP * (float(p["tank_end_C"]) - T_START) / 3.6e6, 2),
                 float(p["electricity_kWh"]), p["tank_end_C"], p["time_45C_UTC"]])   # add PV/T row
    rows.append([day, "PV only", 0.0, round(run_pv(hours, G, Ta), 2), "-", "-"])   # plain PV: electricity only, no heat
    heat, T_end, t45 = run_thermal(hours, G, Ta)             # run the plain thermal collector
    rows.append([day, "Thermal only", round(heat, 2), 0.0, round(T_end, 1), clock(t45)])   # plain thermal: heat only

# Print the results table on screen
print(f"{'Day':<9}{'System':<14}{'Heat kWh':>9}{'Elec kWh':>9}{'Total':>7}{'Tank end C':>11}{'45 C at (UTC)':>15}") # column headings
for r in rows:                                               # one line per system per day
    print(f"{r[0]:<9}{r[1]:<14}{r[2]:>9.2f}{r[3]:>9.2f}{r[2] + r[3]:>7.2f}{str(r[4]):>11}{r[5]:>15}")   # total = heat + electricity

# Save the same table as a CSV file
with open("comparison_results.csv", "w") as f:               # create or overwrite the file
    f.write("day,system,heat_kWh,electricity_kWh,total_kWh,tank_end_C,time_45C_UTC\n")   # header row
    for r in rows:                                           # one line per system per day
        f.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{round(r[2] + r[3], 2)},{r[4]},{r[5]}\n")

# Bar chart: heat and electricity for each system on each day
fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for ax, day in zip(axes, DAYS):
    sel = [r for r in rows if r[0] == day]
    names = [r[1] for r in sel]
    ax.bar(names, [r[3] for r in sel], label="Electricity")
    ax.bar(names, [r[2] for r in sel], bottom=[r[3] for r in sel], label="Heat to tank")
    ax.set_title(f"{day} September day")
    ax.grid(axis="y", alpha=0.3)
axes[0].set_ylabel("Energy from 8 m2 (kWh/day)")
axes[0].legend()
plt.tight_layout()
plt.savefig("comparison_plot.png", dpi=200)
print("\nSaved comparison_results.csv and comparison_plot.png")

# Keep the window open when the script is double-clicked.
input("\nPress Enter to close...")
