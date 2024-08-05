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
curr_geometry.rotor.n_channels = 30
curr_options.rotor.n_rotor = 4000

tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStatorMil, rotor=SPRotor)

# Input Constraints Data
P_in = 8000000        # [Pa]
# P_out = 3800000         # [Pa]
T_in_c = 94             # [°C]
T_in = T_in_c + 273.15  # [K]

# From Ravinath Preliminary Analysis
tt.rotor.gap_losses_control = False

# tt.rotor.rpm = 15000
tt.rotor.dv_perc = 0

Z_stat = 1
tt.geometry.stator.Z_stat = Z_stat

throat_width = 0.0035

tt.geometry.throat_width = throat_width
tt.geometry.disc_thickness = 0.0001
tt.geometry.rotor.n_discs = 1
tt.geometry.H_s = 0.00005
tt.geometry.d_main = 0.15
tt.geometry.rotor.b_channel = 0.00005         # [m]
tt.geometry.rotor.d_ratio = 3            # [m]
tt.n_packs = 30
tt.geometry.n_channels = tt.n_packs

# nozzle_width = tt.geometry.disc_thickness * (tt.n_packs - 1) + tt.geometry.rotor.b_channel * tt.n_packs
nozzle_width = tt.geometry.H_s * tt.n_packs
a_width = throat_width * nozzle_width

inlet_thermopoint = tt.points[0].duplicate()
inlet_thermopoint.set_variable("P", P_in)
inlet_thermopoint.set_variable("T", T_in)
inlet_ss = inlet_thermopoint.get_variable("c")

mfr_max = inlet_ss * a_width * inlet_thermopoint.get_variable("rho") * Z_stat

mfr = 0.1                       # [kg/s]
tt.stator.m_dot_s = mfr

inlet_tp_it = tt.points[0].duplicate()
inlet_tp_it.set_variable("p", P_in)
inlet_tp_it.set_variable("h", inlet_thermopoint.get_variable("h"))

rho_it = inlet_tp_it.get_variable("rho")
T_it = inlet_tp_it.get_variable("T")

A = np.pi * tt.geometry.d_main * tt.geometry.rotor.b_channel
ur = mfr / (rho_it * A)

v_nozzle = mfr / (rho_it * a_width * Z_stat)

alpha_in = 85           # [°]
alpha_rad_in = alpha_in * np.pi / 180

# %%------------   CALCULATIONS                            ----------------------------------------------------------> #

tt.fix_rotor_inlet_condition(P_in, T_in, alpha_rad_in, v_nozzle)
rotor_array = tt.rotor.get_rotor_array()
tt.evaluate_rotor_performances()

P_out = rotor_array[-1, 11]

tang_speed = tt.rotor.omega * tt.geometry.d_main / 2

# %%------------             TRAJECTORY PLOT                        -------------------------------------------------> #
fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

ax.plot(rotor_array[:, 24] * np.pi / 180, rotor_array[:, 1])
ax.set_rmax(0.075)
ax.set_rticks([tt.geometry.d_main / tt.geometry.rotor.d_ratio /2,
               (tt.geometry.d_main - tt.geometry.d_main / tt.geometry.rotor.d_ratio) / 2])
ax.grid(True)

ax.set_title("theta_rel")

plt.show()

# %%------------             VELOCITY PLOT                        ---------------------------------------------------> #

plt.rcParams['figure.constrained_layout.use'] = True
fig1, (ax3, ax4) = plt.subplots(1, 2, figsize=(10, 6))

# Absolute Velocities
ax3.plot(rotor_array[:, 1], rotor_array[:, 3], label="Tangential Speed vt")
ax3.plot(rotor_array[:, 1], rotor_array[:, 4], label="Radial Speed vr")
ax3.plot(rotor_array[:, 1], rotor_array[:, 2], label="Rotational Speed u")

ax3.set_xlabel("Radial Distance [m]", color = "Black", fontsize = 14)
ax3.set_ylabel("Velocity [m/s]", color = "Black", fontsize = 14)
ax3.legend(edgecolor = "Black", fontsize = 12)
ax3.grid()

# Relative Velocities
ax4.plot(rotor_array[:, 1], rotor_array[:, 5], label="Tangential Speed wt")
ax4.plot(rotor_array[:, 1], rotor_array[:, 6], label="Radial Speed wr")
ax4.plot(rotor_array[:, 1], rotor_array[:, 2], label="Rotational Speed u")

ax4.set_xlabel("Radial Distance [m]", color = "Black", fontsize = 14)
ax4.set_ylabel("Velocity [m/s]", color = "Black", fontsize = 14)
ax4.legend(edgecolor = "Black", fontsize = 12)
ax4.grid()

plt.show()