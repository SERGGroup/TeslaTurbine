# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.multi_phase import TPRotor, TPTeslaGeometry, TPTeslaOptions, TPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

curr_geometry = TPTeslaGeometry()
curr_geometry.d_main = 0.3                  # [m]
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
# P_out = 427308          # [Pa]
P_out = np.linspace(995000, 50000, 100)

# Turbine Initialization
tt = BaseTeslaTurbine('R1234ze', curr_geometry, curr_options, stator=TPStatorMil, rotor=TPRotor)

tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("x", x_in)
tt.static_points[1].set_variable("P", P_out)

tt.stator.stator_eff = 0.9
tt.options.stator.metastability_check = True

tt.P_in = P_in

output_array = np.empty((len(P_out), 3))

# %%------------           SOLVING                         ----------------------------------------------------------> #
for i in tqdm(range(len(P_out))):

    tt.static_points[1].set_variable("P", P_out[i])
    tt.stator.solve()

    output_array[i, 0] = tt.stator.m_dot_s
    output_array[i, 1] = tt.stator.out_speed
    output_array[i, 2] = tt.stator.x_out

# %%------------           PLOT RESULTS                    ----------------------------------------------------------> #

fig