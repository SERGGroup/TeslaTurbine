# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.single_phase import SPRotor, SPTeslaGeometry, SPTeslaOptions, SPStatorMil, SPStator
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
import matplotlib.pyplot as plt

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

# Input Constraints Data
P_in = 10100000  # [Pa]
P_out = 3900000         # [Pa]
T_in_c = 94  # [°C]
T_in = T_in_c + 273.15  # [K]

m_rif = 0.157         # [kg/s]


# %%------------   CALCULATION                             ----------------------------------------------------------> #

# INITIALIZING TURBINE CONDITIONS
curr_geometry = SPTeslaGeometry()
curr_options = SPTeslaOptions()
curr_options.rotor.integr_variable = 0.04  # [-]
curr_options.stator.metastability_check = False
curr_options.rotor.sp_check = True
curr_options.rotor.n_rotor = 4000

# Main design Parameters
curr_geometry.rotor.b_channel = 0.00005
curr_geometry.rotor.d_ratio = 3  # [m]
curr_geometry.d_main = 0.15  # [m]

tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStatorMil, rotor=SPRotor)

tt.rotor.rpm = 15000

Z_stat = 1
tt.geometry.stator.Z_stat = Z_stat

throat_width = 0.00114
tt.geometry.throat_width = throat_width

tt.geometry.disc_thickness = 0.0001
tt.geometry.rotor.n_discs = 1
tt.geometry.H_s = 0.00015
tt.m_dot_tot = m_rif
tt.geometry.n_channels = tt.n_packs

alpha_in = 89  # [°]
tt.geometry.alpha1 = alpha_in

tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)
tt.stator.stator_eff = 0.9
tt.rotor.gap_losses_control = True

tt.P_in = P_in
tt.P_out = P_out
tt.T_in = T_in

tt.iterate_pressure()
tt.evaluate_performances()
rotor_array = tt.rotor.get_rotor_array()


# %%------------             TRAJECTORY PLOT                        -------------------------------------------------> #
fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={"projection": "polar"})

ax1.plot(rotor_array[:, 24] * np.pi / 180, rotor_array[:, 1], zorder=1, color='Darkblue', linewidth='1.3')
ax1.set_rmax(0.075)
ax1.set_rticks([tt.geometry.d_main / tt.geometry.rotor.d_ratio /2,
               (tt.geometry.d_main - tt.geometry.d_main / tt.geometry.rotor.d_ratio) / 2])
ax1.tick_params(axis='y', labelsize=9, zorder=10)
ax1.tick_params(axis='x', labelsize=9, zorder=10)
ax1.grid(True)

ax1.set_title("Absolute Reference Frame")

ax2.plot(rotor_array[:, 25] * np.pi / 180, rotor_array[:, 1], color='Darkred', linewidth='1.3')
ax2.set_rmax(0.075)
ax2.set_rticks([tt.geometry.d_main / tt.geometry.rotor.d_ratio /2,
               (tt.geometry.d_main - tt.geometry.d_main / tt.geometry.rotor.d_ratio) / 2])
ax2.tick_params(axis='y', labelsize=9)
ax2.tick_params(axis='x', labelsize=9)
ax2.grid(True)

ax2.set_title("Relative Reference Frame")

fig.tight_layout()
plt.show()


# %%------------             VELOCITY PLOT                        ---------------------------------------------------> #

plt.rcParams['figure.constrained_layout.use'] = True
fig1, (ax3, ax4) = plt.subplots(1, 2, figsize=(9, 5))

# Absolute Velocities
ax3.plot(rotor_array[1:, 1], rotor_array[1:, 3], label="Tangential Speed Vt", color = "Darkred")
ax3.plot(rotor_array[1:, 1], rotor_array[1:, 4], label="Radial Speed Vr", color = "Darkgreen")
ax3.plot(rotor_array[1:, 1], rotor_array[1:, 2], label="Rotational Speed U", color = "Darkblue")

ax3.set_xlim(0.025, 0.075)
ax3.set_ylim(0, 225)

ax3.set_xlabel("Radial Position [m]", color="Black", fontsize=14)
ax3.set_ylabel("Velocity [m/s]", color="Black", fontsize=14)
ax3.legend(edgecolor="Black", fontsize=12)
ax3.grid()

# Pressure
ax4.plot(rotor_array[1:, 1], rotor_array[1:, 11] / 100000, label="Static Pressure", color="black")

ax4.set_xlabel("Radial Position [m]", color = "Black", fontsize = 14)
ax4.set_ylabel("Static Pressure [bar]", color = "Black", fontsize = 14)

ax4.set_xlim(0.025, 0.075)
ax4.set_ylim(38, 50)
ax4.grid()

plt.show()

# %%------------             VELOCITY PLOT                        ---------------------------------------------------> #
fig2, ax = plt.subplots(figsize=(5,4))

ax.plot(rotor_array[1:, 1], rotor_array[1:, 3], label="Tangential Speed Vt", color = "Darkred")
ax.plot(rotor_array[1:, 1], rotor_array[1:, 4], label="Radial Speed Vr", color = "Darkgreen")
ax.plot(rotor_array[1:, 1], rotor_array[1:, 2], label="Rotational Speed U", color = "Darkblue")

ax.set_xlim(0.025, 0.075)
ax.set_ylim(0, 225)

ax.set_xlabel("Radial Position [m]", color="Black", fontsize=14)
ax.set_ylabel("Velocity [m/s]", color="Black", fontsize=14)
ax.legend(edgecolor="Black", fontsize=12)
ax.grid()

plt.show()