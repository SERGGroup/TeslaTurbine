# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
import numpy as np
from main_code.sub_classes.single_phase import SPRotor, SPTeslaGeometry, SPTeslaOptions, SPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import matplotlib.pyplot as plt
import pandas as pd

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

curr_geometry = SPTeslaGeometry()
curr_options = SPTeslaOptions()
curr_options.rotor.integr_variable = 0.04           # [-]
curr_options.stator.metastability_check = False
curr_options.rotor.sp_check = True
width = 0.02                                        # [m]
disc_thickness = 0.0001                             # [m]
curr_options.rotor.n_rotor = 3000

# Main design Parameters
curr_geometry.rotor.b_channel = 0.0003         # [m]
curr_geometry.rotor.d_ratio = 3.75            # [m]
curr_geometry.d_main = 0.3                    # [m]

n_channels = (width + disc_thickness) / (curr_geometry.rotor.b_channel + disc_thickness)
curr_geometry.rotor.n_channels = n_channels

tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStatorMil, rotor=SPRotor)

# Input Constraints Data
P_in = 10000000         # [Pa]
P_out = 3800000         # [Pa]
T_in_c = 94             # [°C]
T_in = T_in_c + 273.15  # [K]

# From Ravinath Preliminary Analysis
tt.rotor.gap_losses_control = False
Z_stat = 4
tt.geometry.stator.Z_stat = Z_stat

alpha_in = 85           # [°]
alpha_rad_in = alpha_in * np.pi / 180
throat_width = 0.0005
nozzle_width = width
a_width = throat_width * nozzle_width


inlet_thermopoint = tt.points[0].duplicate()
inlet_thermopoint.set_variable("P", P_in)
inlet_thermopoint.set_variable("T", T_in)
inlet_ss = inlet_thermopoint.get_variable("c")

mfr_max = inlet_ss * a_width * inlet_thermopoint.get_variable("rho") * Z_stat

dv_perc_array = np.linspace(-0.5, 0.5, 20)
results_array = np.empty((len(dv_perc_array), 5))

# %%------------   CALCULATIONS                            ----------------------------------------------------------> #
it = 20

for j in range(len(dv_perc_array)):

    v_nozzle = 0.
    tt.rotor.dv_perc = dv_perc_array[j]

    mfr_up = mfr_max
    mfr_down = 0.

    print(' ')
    print('Starting {}-th Iteration on  Mass Flow Rate -- Limit Values are: MIN {} , MAX {}'.format(j, 0, mfr_max))

    for i in range(it):

        mfr = (mfr_up + mfr_down) / 2
        tt.stator.m_dot_s = mfr
        v_nozzle = mfr / (a_width * inlet_thermopoint.get_variable("rho") * Z_stat)

        tt.fix_rotor_inlet_condition(P_in, T_in, alpha_rad_in, v_nozzle)
        rotor_array = tt.rotor.get_rotor_array()
        tt.evaluate_rotor_performances()

        outlet_pressure = rotor_array[-1, 11]
        error = abs((outlet_pressure - P_out) / P_out)

        if error < 0.0001:
            break
        elif P_out < outlet_pressure:
            mfr_down = mfr
        else:
            mfr_up = mfr

        if abs(mfr - mfr_max) < 0.01:

            print('WARNING!!! The NOZZLE is CHOKED!')
            print(' ')
            break

        print(mfr)

    results_array[j, 0] = tt.Eta_rotor_ss
    results_array[j, 1] = tt.power
    results_array[j, 2] = v_nozzle / inlet_ss
    results_array[j, 3] = tt.rotor.rpm
    results_array[j, 4] = tt.power * tt.geometry.rotor.n_channels

n = 0
for i in range(len(dv_perc_array) - n):

    if np.isnan(results_array[i, 0]) or results_array[i, 1] < 0:

        results_array = np.delete(results_array, i, 0)
        n += 1

# df = pd.DataFrame(data=results_array, columns=['Eta [-], Power [W/ch], Stator Outlet Mach Number [-], RPM [rpm], Power [W]'])

# %%------------             TRAJECTORY PLOT                        -------------------------------------------------> #
fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

ax.plot(rotor_array[:, 25] * np.pi / 180, rotor_array[:, 1])
ax.set_rmax(0.15)
ax.set_rticks([0.025, 0.05, 0.075, 0.1, 0.125])
ax.set_rlabel_position(-22.5)
ax.grid(True)

ax.set_title("theta_rel")

plt.show()

# %%------------             RPM VARIATION PLOT                     -------------------------------------------------> #

fig1, ax1 = plt.subplots()

ax1.set_xlabel('RPM [rpm]', fontsize=14)
ax1.set_ylabel('Power [W]', fontsize=14)
ax1.plot(results_array[:, 3], results_array[:, 4], color='Darkblue')

ax2 = ax1.twinx()

ax2.set_ylabel('Efficiency [-]', fontsize=14)
ax2.plot(results_array[:, 3], results_array[:, 0], color='Darkred')

fig1.tight_layout()
plt.show()