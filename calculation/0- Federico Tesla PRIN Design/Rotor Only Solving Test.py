# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
import numpy as np
from main_code.sub_classes.single_phase import SPRotor, SPTeslaGeometry, SPTeslaOptions, SPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import matplotlib.pyplot as plt

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

curr_geometry = SPTeslaGeometry()
curr_options = SPTeslaOptions()
curr_options.rotor.integr_variable = 0.04           # [-]
curr_options.stator.metastability_check = False
curr_options.rotor.sp_check = True
curr_geometry.rotor.n_channels = 100
curr_options.rotor.n_rotor = 3000

# Main design Parameters
curr_geometry.rotor.b_channel = 0.0001         # [m]
curr_geometry.rotor.d_ratio = 3.75            # [m]
curr_geometry.d_main = 0.3                    # [m]

tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStatorMil, rotor=SPRotor)

# Input Constraints Data
P_in = 10000000         # [Pa]
P_out = 3800000         # [Pa]
T_in_c = 94             # [°C]
T_in = T_in_c + 273.15  # [K]

# From Ravinath Preliminary Analysis
tt.rotor.gap_losses_control = False
tt.rotor.rpm = 7500
Z_stat = 4
tt.geometry.stator.Z_stat = Z_stat

alpha_in = 87           # [°]
alpha_rad_in = alpha_in * np.pi / 180
throat_width = 0.0005
nozzle_width = 0.02
a_width = throat_width * nozzle_width

inlet_thermopoint = tt.points[0].duplicate()
inlet_thermopoint.set_variable("P", P_in)
inlet_thermopoint.set_variable("T", T_in)
inlet_ss = inlet_thermopoint.get_variable("c")

mfr_max = inlet_ss * a_width * inlet_thermopoint.get_variable("rho") * Z_stat

# %%------------   CALCULATIONS                            ----------------------------------------------------------> #
mfr_up = mfr_max
mfr_down = 0.
it = 100

print(' ')
print('Iteration on Mass Flow Rate -- Limit Values are: MIN {} , MAX {}'.format(0, mfr_max))

for i in range(it):

    mfr = (mfr_up + mfr_down) / 2
    tt.stator.m_dot_s = mfr
    v_nozzle = mfr / (a_width * inlet_thermopoint.get_variable("rho") * Z_stat)

    tt.fix_rotor_inlet_condition(P_in, T_in, alpha_rad_in, v_nozzle)
    rotor_array = tt.rotor.get_rotor_array()
    tt.evaluate_rotor_performances()

    outlet_pressure = rotor_array[-1, 11]
    error = abs((outlet_pressure - P_out) / P_out)
    abs_error = abs(outlet_pressure - P_out)

    if error < 0.0001:

        break

    elif P_out < outlet_pressure:

        mfr_down = mfr

    else:

        mfr_up = mfr

    # rpm = tt.rotor.rpm
    dv_perc = tt.rotor.dv_perc

    if abs(mfr - mfr_max) < 0.001:

        print('WARNING!!! The NOZZLE is CHOKED')
        print(' ')

        break

    print('Mass Flow Rate is {}'.format(mfr))
    # print('rpm is {}'.format(rpm))
    print('Velocity Ratio is {}'.format(dv_perc))
    print(' ')


# %%------------             TRAJECTORY PLOT                        -------------------------------------------------> #
fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

ax.plot(rotor_array[:, 24] * np.pi / 180, rotor_array[:, 1])
ax.set_rmax(0.15)
ax.set_rticks([0.025, 0.05, 0.075, 0.1, 0.125])
ax.set_rlabel_position(-22.5)
ax.grid(True)

ax.set_title("theta_rel")

plt.show()

# %%------------             VELOCITY PLOT                        ---------------------------------------------------> #

plt.rcParams['figure.constrained_layout.use'] = True
fig1, (ax3, ax4) = plt.subplots(1, 2, figsize=(10, 6))

# Absolute Velocities
ax3.plot(rotor_array[:, 1], rotor_array[:, 3], label="Tangential Speed vt")
ax3.plot(rotor_array[:, 1], rotor_array[:, 4], label="Radial Speed vr")
ax3.plot(rotor_array[:, 1], rotor_array[:, 2], label="Rotative Speed u")

ax3.set_xlabel("Radial Distance [m]", color = "Black", fontsize = 14)
ax3.set_ylabel("Velocity [m/s]", color = "Black", fontsize = 14)
# ax3.set_xlim(0.04, 0.1)
# ax3.set_ylim(0, 75)
ax3.legend(edgecolor = "Black", fontsize = 12)
ax3.grid()

# Relative Velocities
ax4.plot(rotor_array[:, 1], rotor_array[:, 5], label="Tangential Speed wt")
ax4.plot(rotor_array[:, 1], rotor_array[:, 6], label="Radial Speed wr")
ax4.plot(rotor_array[:, 1], rotor_array[:, 2], label="Rotative Speed u")

ax4.set_xlabel("Radial Distance [m]", color = "Black", fontsize = 14)
ax4.set_ylabel("Velocity [m/s]", color = "Black", fontsize = 14)
# ax4.set_xlim(0.04, 0.1)
# ax4.set_ylim(0, 25)
ax4.legend(edgecolor = "Black", fontsize = 12)
ax4.grid()

plt.show()
