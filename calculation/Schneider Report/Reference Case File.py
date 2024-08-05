# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.multi_phase import TPRotor, TPTeslaGeometry, TPTeslaOptions, TPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

curr_geometry = TPTeslaGeometry()
curr_geometry.d_main = 0.2
curr_geometry.rotor.b_channel = 0.0005  # [m]
curr_geometry.disc_thickness = 0.0008  # [m]
curr_geometry.alpha1 = 85  # [°]
curr_geometry.stator.Z_stat = 4  # [-]

curr_geometry.rotor.gap = 0.001  # [m]

curr_geometry.throat_width = 0.015              # [m]
curr_geometry.rotor.n_discs = 1                 # [-]
curr_geometry.rotor.d_ratio = 4                 # [-]
curr_geometry.rotor.n_packs = 1                 # [-]
curr_geometry.rotor.roughness = 0.0000005       # [m]

single_hs = (curr_geometry.rotor.n_discs - 1) * curr_geometry.disc_thickness + curr_geometry.rotor.n_discs * curr_geometry.rotor.b_channel
curr_geometry.H_s = single_hs * curr_geometry.rotor.n_packs

# Solving Parameters
curr_options = TPTeslaOptions()
curr_options.rotor.integr_variable = 0.03
curr_options.rotor.n_rotor = 400

tt = BaseTeslaTurbine('R1234ze', curr_geometry, curr_options, stator=TPStatorMil, rotor=TPRotor)

n_setpoints = 10

x_in = 0                                        # [-]
P_out = 427308                                  # [Pa]
P_in = 997233

rpm_opt = 4224

tt.rotor.rpm = rpm_opt
tt.points[0].set_variable("p", P_in)
tt.points[0].set_variable("x", x_in)
tt.n_packs = 198
tt.stator.stator_eff = 0.9
tt.rotor.gap_losses_control = True

tt.options.stator.metastability_check = True
tt.options.rotor.profile_rotor = True
tt.options.rotor.sp_check = False
tt.options.rotor.tp_epsilon_model = "sarti"

tt.P_in = P_in
tt.P_out = P_out

# %%------------   CALCULATIONS                            ----------------------------------------------------------> #

tt.iterate_pressure()
tt.evaluate_performances()
rotor_array = tt.rotor.get_rotor_array()

# %%------------             TRAJECTORY PLOT                        -------------------------------------------------> #
fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

ax.plot(rotor_array[:, 25] * np.pi / 180, rotor_array[:, 1])
ax.set_rmax(0.1)
ax.set_rticks([0.025, 0.05, 0.075, 0.1])
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