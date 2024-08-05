# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
import numpy as np
from main_code.sub_classes.single_phase import SPRotor, SPTeslaGeometry, SPTeslaOptions, SPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import matplotlib.pyplot as plt
from REFPROPConnector import ThermodynamicPoint as TP

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

curr_options = SPTeslaOptions()
curr_options.rotor.integr_variable = 0.04           # [-]
curr_options.stator.metastability_check = False
curr_options.rotor.sp_check = True
curr_options.rotor.n_rotor = 3000

# Input Constraints Data
P_in = 10000000         # [bar]
P_out = 3800000         # [bar]
T_in_c = 94             # [°C]
T_in = T_in_c + 273.15  # [K]

# From Ravinath Preliminary Analysis
alpha_in = 75           # [°]
alpha_rad_in = alpha_in * np.pi / 180
throat_width = 0.0005
nozzle_width = 0.02
a_width = throat_width * nozzle_width
Z_stat = 4

inlet_thermopoint = TP(['CarbonDioxide'], [1], unit_system="MASS BASE SI")
inlet_thermopoint.set_variable("P", P_in)
inlet_thermopoint.set_variable("T", T_in)
inlet_ss = inlet_thermopoint.get_variable("c")

mfr_max = inlet_ss * a_width * inlet_thermopoint.get_variable("rho") * Z_stat
# mfr_max = 12

b_channel_array = np.linspace(0.0001, 0.0009, 10)
Eta = np.empty(len(b_channel_array))
Power = np.empty(len(b_channel_array))
Mach_number = np.empty(len(b_channel_array))

output_array = np.empty((len(b_channel_array), 8))
output_array[:] = np.NaN

Sound_speed = inlet_thermopoint.get_variable("c")
Mach = 0.7

# %%------------   CALCULATIONS                            ----------------------------------------------------------> #
for j in range(len(b_channel_array)):

    v_nozzle = 0.
    curr_geometry = SPTeslaGeometry()
    curr_geometry.rotor.b_channel = b_channel_array[j]
    curr_geometry.rotor.d_ratio = 3.75                   # [m]
    curr_geometry.d_main = 0.3                           # [m]
    curr_geometry.stator.Z_stat = Z_stat

    tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStatorMil, rotor=SPRotor)

    tt.rotor.gap_losses_control = False
    tt.rotor.rpm = 7000
    # tt.rotor.dv_perc = 0
    tt.geometry.rotor.n_channels = nozzle_width / tt.geometry.rotor.b_channel

    mfr_up = mfr_max
    mfr_down = 0.
    it = 20
    i = 0.

    print(' ')
    print('Starting {}-th Iteration on  Mass Flow Rate -- Limit Values are: MIN {} , MAX {}'.format(j, 0, mfr_max))

    for i in range(it):

        mfr = (mfr_up + mfr_down) / 2
        tt.stator.m_dot_s = mfr
        # v_nozzle = mfr / (a_width * inlet_thermopoint.get_variable("rho") * Z_stat)
        v_nozzle = Mach * Sound_speed

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

        if i == it:

            print('WARNING!!! The SIMULATION is NOT CONVERGING!')
            print(' ')

    if i == it or tt.Eta_rotor_ss > 1 or tt.Eta_rotor_ss < 0:

        pass

    else:

        output_array[j, 0] = tt.Eta_rotor_ss
        output_array[j, 1] = tt.work
        output_array[j, 2] = v_nozzle / inlet_ss
        output_array[j, 3] = tt.rotor.rotor_points[-1].speed.v
        output_array[j, 4] = tt.power
        output_array[j, 5] = tt.geometry.rotor.n_channels
        output_array[j, 6] = tt.rotor.m_dot_ch
        output_array[j, 7] = tt.power * tt.geometry.rotor.n_channels

    # Eta[j] = tt.Eta_rotor_ss
    # Power[j] = tt.power
    # Mach_number[j] = v_nozzle / inlet_ss

# %%------------             TRAJECTORY PLOT                        -------------------------------------------------> #
fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

ax.plot(rotor_array[:, 24] * np.pi / 180, rotor_array[:, 1])
ax.set_rmax(0.15)
ax.set_rticks([0.05, 0.1])
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
ax3.legend(edgecolor = "Black", fontsize = 12)
ax3.grid()

# Relative Velocities
ax4.plot(rotor_array[:, 1], rotor_array[:, 5], label="Tangential Speed wt")
ax4.plot(rotor_array[:, 1], rotor_array[:, 6], label="Radial Speed wr")
ax4.plot(rotor_array[:, 1], rotor_array[:, 2], label="Rotative Speed u")

ax4.set_xlabel("Radial Distance [m]", color = "Black", fontsize = 14)
ax4.set_ylabel("Velocity [m/s]", color = "Black", fontsize = 14)
ax4.legend(edgecolor = "Black", fontsize = 12)
ax4.grid()

plt.show()
