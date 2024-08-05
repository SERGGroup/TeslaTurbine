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
curr_options.rotor.n_rotor = 4000

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

RPM = np.linspace(5000,10000, 11)
# tt.rotor.rpm = 7000
# tt.rotor.dv_perc = 0.

Z_stat = 4
tt.geometry.stator.Z_stat = Z_stat

# alpha_in = 87           # [°]
# alpha_rad_in = alpha_in * np.pi / 180
throat_width = 0.0005
nozzle_width = 0.02
a_width = throat_width * nozzle_width

tt.geometry.throat_width = throat_width
tt.geometry.disc_thickness = 0.0001
tt.geometry.rotor.n_discs = 1
tt.geometry.H_s = 0.0001
tt.geometry.rotor.n_packs = tt.geometry.rotor.n_discs

inlet_thermopoint = tt.points[0].duplicate()
inlet_thermopoint.set_variable("P", P_in)
inlet_thermopoint.set_variable("T", T_in)
inlet_ss = inlet_thermopoint.get_variable("c")

mfr_max = inlet_ss * a_width * inlet_thermopoint.get_variable("rho") * Z_stat

mfr = 0.1                           # [kg/s]
tt.stator.m_dot_s = mfr

output_array = np.empty((len(RPM), 7))
output_array[:] = np.NaN

# %%------------   CALCULATIONS                            ----------------------------------------------------------> #

print(' ')
print('Iteration on Inlet Pressure -- Limit Values are: MIN {} , MAX {}'.format(P_in, P_out))
for j in range(len(RPM)):

    print(' ')
    print('Iteration on {}-th RPM  -- RPM = {}'.format(j+1, RPM[j]))
    print(' ')

    P_up = P_in
    P_down = P_out
    it = 30
    P_it = 0.

    for i in range(it):

        P_it = (P_up + P_down) / 2

        print('Starting {}-th iteration on pressure -- P_it = {}'.format(i+1, P_it))

        inlet_tp_it = tt.points[0].duplicate()
        inlet_tp_it.set_variable("p", P_it)
        inlet_tp_it.set_variable("h", inlet_thermopoint.get_variable("h"))

        rho_it = inlet_tp_it.get_variable("rho")
        T_it = inlet_tp_it.get_variable("T")

        A = np.pi * tt.geometry.d_main * tt.geometry.rotor.b_channel
        ur = mfr / (rho_it * A)

        # v_nozzle = 0.8 * inlet_tp_it.get_variable("c")
        v_nozzle = mfr / (rho_it * a_width)

        tt.rotor.rpm = RPM[j]
        tt.fix_rotor_inlet_condition2(P_it, T_it, ur, v_nozzle)
        rotor_array = tt.rotor.get_rotor_array()
        tt.evaluate_rotor_performances()

        outlet_pressure = rotor_array[-1, 11]
        error = abs((outlet_pressure - P_out) / P_out)
        abs_error = abs(outlet_pressure - P_out)

        if error < 0.0001:

            break

        elif P_out > outlet_pressure:

            P_down = P_it

        else:

            P_up = P_it

    output_array[j, 0] = tt.Eta_rotor_ss
    output_array[j, 1] = tt.work
    output_array[j, 2] = P_it
    output_array[j, 3] = tt.rotor.rotor_points[-1].speed.v
    output_array[j, 4] = tt.power
    output_array[j, 5] = tt.rotor.rotor_inlet_speed.alpha
    output_array[j, 6] = tt.work2

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