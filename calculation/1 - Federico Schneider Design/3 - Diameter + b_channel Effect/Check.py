# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.multi_phase import TPRotor, TPTeslaGeometry, TPTeslaOptions, TPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

curr_geometry = TPTeslaGeometry()
curr_geometry.d_main = 0.295167               # [m]
curr_geometry.throat_width = 0.003          # [m]
curr_geometry.disc_thickness = 0.0008       # [m]
curr_geometry.alpha1 = 85                   # [°]
curr_geometry.stator.Z_stat = 4             # [-]

curr_geometry.rotor.gap = 0.001             # [m]

curr_geometry.rotor.b_channel = 0.0001      # [m]
curr_geometry.rotor.n_discs = 1             # [-]
curr_geometry.rotor.d_ratio = 3.5           # [-]
curr_geometry.rotor.n_packs = 1             # [-]
curr_geometry.rotor.roughness = 0.0000005   # [m]

single_hs = ( curr_geometry.rotor.n_discs - 1) * curr_geometry.disc_thickness + curr_geometry.rotor.n_discs * curr_geometry.rotor.b_channel
curr_geometry.H_s = single_hs * curr_geometry.rotor.n_packs

# Solving Parameters
curr_options = TPTeslaOptions()
curr_options.rotor.integr_variable = 0.03
curr_options.rotor.n_rotor = 1200

# Thermodynamic Parameters
P_in = 997233           # [Pa]
x_in = 0                # [-]
P_out = 427308          # [Pa]
m_refr = 3.336  # [kg/s]

# Turbine Initialization
tt = BaseTeslaTurbine('R1234ze', curr_geometry, curr_options, stator=TPStatorMil, rotor=TPRotor)

tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("x", x_in)
tt.stator.stator_eff = 0.9
tt.rotor.gap_losses_control = True

tt.rotor.dv_perc = 0
tt.options.stator.metastability_check = True

tt.P_in = P_in
# tt.static_points[1].set_variable("p", 658000)
tt.P_out = P_out

# %%------------   CALCULATION                             ----------------------------------------------------------> #

# tt.stator.solve()
# tt.rotor.solve()
tt.iterate_pressure()
rotor_array = tt.rotor.get_rotor_array()
tt.evaluate_performances()

tt.geometry.rotor.n_channels = np.rint(m_refr / tt.stator.m_dot_s)

Power = tt.power * tt.geometry.rotor.n_channels

# %%------------             TRAJECTORY PLOT                        -------------------------------------------------> #
fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

ax.plot(rotor_array[:, 24] * np.pi / 180, rotor_array[:, 1])
ax.set_rmax(0.09625)
ax.set_rticks([0.025, 0.05, 0.075])
ax.set_rlabel_position(-22.5)
ax.grid(True)

ax.set_title("theta_rel")

plt.show()