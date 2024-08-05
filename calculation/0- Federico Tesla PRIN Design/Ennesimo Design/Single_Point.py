# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.single_phase import SPRotor, SPTeslaGeometry, SPTeslaOptions, SPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
import matplotlib.pyplot as plt

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

# Input Constraints Data
P_in = 10000000  # [Pa]
P_out = 3800000         # [Pa]
T_in_c = 94  # [°C]
T_in = T_in_c + 273.15  # [K]

m_rif = 0.1             # [kg/s]


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

tt.rotor.gap_losses_control = False
tt.rotor.rpm = 15000

Z_stat = 1
tt.geometry.stator.Z_stat = Z_stat

throat_width = 0.001
tt.geometry.throat_width = throat_width

tt.geometry.disc_thickness = 0.0001
tt.geometry.rotor.n_discs = 1
tt.geometry.H_s = 0.00005
tt.m_dot_tot = m_rif
tt.geometry.n_channels = tt.n_packs

alpha_in = 85  # [°]
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
fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={"projection": "polar"}, figsize=(10, 6))

ax1.plot(rotor_array[:, 24] * np.pi / 180, rotor_array[:, 1])
ax1.set_rmax(0.075)
ax1.set_rticks([tt.geometry.d_main / tt.geometry.rotor.d_ratio /2,
               (tt.geometry.d_main - tt.geometry.d_main / tt.geometry.rotor.d_ratio) / 2])
ax1.grid(True)

ax1.set_title("Absolute Reference Frame")

ax2.plot(rotor_array[:, 25] * np.pi / 180, rotor_array[:, 1])
ax2.set_rmax(0.075)
ax2.set_rticks([tt.geometry.d_main / tt.geometry.rotor.d_ratio /2,
               (tt.geometry.d_main - tt.geometry.d_main / tt.geometry.rotor.d_ratio) / 2])
ax2.grid(True)

ax2.set_title("Relative Reference Frame")

plt.show()

# %%------------             VELOCITY PLOT                        ---------------------------------------------------> #

plt.rcParams['figure.constrained_layout.use'] = True
fig1, (ax3, ax4) = plt.subplots(1, 2, figsize=(10, 6))

# Absolute Velocities
ax3.plot(rotor_array[:, 1], rotor_array[:, 3], label="Tangential Speed vt", color = "Darkred")
ax3.plot(rotor_array[:, 1], rotor_array[:, 4], label="Radial Speed vr", color = "Darkgreen")
ax3.plot(rotor_array[:, 1], rotor_array[:, 2], label="Rotational Speed u", color = "Darkblue")

ax3.set_xlabel("Radial Distance [m]", color = "Black", fontsize = 14)
ax3.set_ylabel("Velocity [m/s]", color = "Black", fontsize = 14)
ax3.legend(edgecolor = "Black", fontsize = 12)
ax3.grid()
ax3.set_title("Absolute Frame")

# Relative Velocities
ax4.plot(rotor_array[:, 1], rotor_array[:, 5], label="Tangential Speed wt", color = "Darkred")
ax4.plot(rotor_array[:, 1], rotor_array[:, 6], label="Radial Speed wr", color = "Darkgreen")
ax4.plot(rotor_array[:, 1], rotor_array[:, 2], label="Rotational Speed u", color = "Darkblue")

ax4.set_xlabel("Radial Distance [m]", color = "Black", fontsize = 14)
ax4.set_ylabel("Velocity [m/s]", color = "Black", fontsize = 14)
ax4.legend(edgecolor = "Black", fontsize = 12)
ax4.grid()
ax4.set_title("Relative Frame")

plt.show()