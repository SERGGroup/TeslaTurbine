# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.single_phase import SPStator, SPRotor, SPTeslaGeometry, SPTeslaOptions, SPStatorMil
from main_code.base_classes import BaseTeslaTurbine, CALCULATION_FOLDER
import matplotlib.pyplot as plt
import numpy as np

# %%------------       INPUT DATA                         -----------------------------------------------------------> #

# Geometric Parameters
curr_geometry = SPTeslaGeometry()
curr_geometry.d_main = 0.217                        # [m]
curr_geometry.throat_width = 0.001                  # [m]
curr_geometry.disc_thickness = 0.0008               # [m]
curr_geometry.alpha1 = 85                           # [°]
curr_geometry.stator.Z_stat = 4                     # [-]

curr_geometry.rotor.gap = 0.001                     # [m]

curr_geometry.rotor.n_discs = 2                     # [-]
curr_geometry.rotor.b_channel = 0.0001              # [m]
curr_geometry.rotor.n_channels = 60                 # [-]
curr_geometry.rotor.d_ratio = 3.92727               # [-]
curr_geometry.rotor.n_packs = 30                    # [-]

single_hs = (curr_geometry.rotor.n_discs - 1) * curr_geometry.disc_thickness + curr_geometry.rotor.n_discs * curr_geometry.rotor.b_channel
curr_geometry.H_s = single_hs * curr_geometry.rotor.n_packs

# Solving Parameters
curr_options = SPTeslaOptions()
curr_options.rotor.integr_variable = 0.03
curr_options.rotor.n_rotor = 400

# Thermodynamic Parameters
P_in = 644701           # [Pa]
T_in = 87.13 + 273.15   # [K]
P_out = 405224          # [Pa]

# # Thermodynamic Parameters
# P_in = 473535           # [Pa]
# T_in = 73.42 + 273.15   # [K]
# P_out = 311123          # [Pa]

#Tesla Turbine Initialization
tt = BaseTeslaTurbine("R1233zd(E)", curr_geometry, curr_options, stator=SPStatorMil, rotor=SPRotor)

tt.P_in = P_in
tt.T_in = T_in
tt.P_out = P_out

tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)
tt.points[3].set_variable("P", P_out)

tt.rotor.rpm = 2500
tt.stator.stator_eff = 0.9
tt.rotor.options.sp_check = True


# %%------------       CALCULATIONS                       -----------------------------------------------------------> #

tt.iterate_total_pressure()
rotor_array = tt.rotor.get_rotor_array()
tt.evaluate_performances()

ss = tt.points[1].get_variable("c")
m_dot = tt.stator.m_dot_s
power = tt.power
rpm = tt.rotor.rpm
torque = power / (rpm * 2 * np.pi / 60)

# %%------------             TRAJECTORY PLOT                        -------------------------------------------------> #
fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

ax.plot(rotor_array[:, 24] * np.pi / 180, rotor_array[:, 1])
ax.set_rmax(0.107)
ax.set_rticks([0.025, 0.05, 0.075])
ax.set_rlabel_position(-22.5)
ax.grid(True)

ax.set_title("theta_rel")

plt.show()

# %%--------------------         PLOT RESULTS             ----------------------------------------------------------> #

fig, (ax1, ax2) = plt.subplots(1, 2, constrained_layout=True)

ax1.plot(rotor_array[:, 1], rotor_array[:, 3], label="vt", c='Darkblue')
ax1.plot(rotor_array[:, 1], rotor_array[:, 2], label="u", c='Darkred')
ax1.plot(rotor_array[:, 1], rotor_array[:, 4], label="vr", c='Darkgreen')
ax1.plot(rotor_array[:, 1], rotor_array[:, 9], label="v", c='black')
ax1.grid()
ax1.legend()
ax1.set_title('Absolute Velocities')
ax1.set_xlabel('Radial Distance [m]')
ax1.set_ylabel('Velocity [m/s]')

ax2.plot(rotor_array[:, 1], rotor_array[:, 5], label="wt", c='Darkblue')
ax2.plot(rotor_array[:, 1], rotor_array[:, 2], label="u", c='Darkred')
ax2.plot(rotor_array[:, 1], rotor_array[:, 6], label="wr", c='Darkgreen')
ax2.plot(rotor_array[:, 1], rotor_array[:, 10], label="w", c='black')
ax2.grid()
ax2.legend()
ax2.set_title('Relative Velocities')
ax2.set_xlabel('Radial Distance [m]')
ax2.set_ylabel('Velocity [m/s]')

plt.show()