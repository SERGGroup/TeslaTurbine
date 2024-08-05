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
P_in = 473535           # [Pa]
T_in = 73.42 + 273.15   # [K]

#Tesla Turbine Initialization
tt = BaseTeslaTurbine("R1233zd(E)", curr_geometry, curr_options, stator=SPStatorMil, rotor=SPRotor)

P_out = np.linspace(473000, 250000, 100)

tt.P_in = P_in
tt.T_in = T_in

tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)

tt.rotor.rpm = 2500
tt.stator.stator_eff = 1
tt.rotor.options.sp_check = True

output_array = np.empty((len(P_out), 3))

#Tesla Turbine Initialization
tt1 = BaseTeslaTurbine("R1233zd(E)", curr_geometry, curr_options, stator=SPStator, rotor=SPRotor)

tt1.P_in = P_in
tt1.T_in = T_in

tt1.points[0].set_variable("P", P_in)
tt1.points[0].set_variable("T", T_in)

tt1.rotor.rpm = 2500
tt1.stator.stator_eff = 0.9
tt1.rotor.options.sp_check = True
output_array1 = np.empty((len(P_out), 3))

# %%------------       CALCULATIONS                       -----------------------------------------------------------> #

for i in range(len(P_out)):
    tt.P_out = P_out[i]
    tt.static_points[1].set_variable("P", P_out[i])

    tt.stator.solve()

    output_array[i, 0] = tt.stator.m_dot_s
    output_array[i, 1] = tt.stator.out_speed
    output_array[i, 2] = tt.points[1].get_variable("c")

    tt1.P_out = P_out[i]
    tt1.static_points[1].set_variable("P", P_out[i])

    tt1.stator.solve()

    output_array1[i, 0] = tt1.stator.m_dot_s
    output_array1[i, 1] = tt1.stator.v_out
    output_array1[i, 2] = tt1.points[1].get_variable("c")

# %%---------------------         PLOT RESULTS (1)             ------------------------------------------------------> #

fig, axs = plt.subplots(1, 2, constrained_layout=True)

axs[0].plot(P_out, output_array[:, 1], label='v_out', c='darkblue')
axs[0].plot(P_out, output_array[:, 2], label='sound_speed', c='darkred')
axs[0].legend()
axs[0].grid()

axs[1].plot(P_out, output_array[:, 0], label='m_dot', c='black')
axs[1].grid()

plt.show()

# %%---------------------         PLOT RESULTS (2)             ------------------------------------------------------> #

fig2, axs2 = plt.subplots(1, 2, constrained_layout=True)

axs2[0].plot(P_out, output_array[:, 1], label='v_out - Mil', c='darkblue')
axs2[0].plot(P_out, output_array[:, 2], label='sound_speed - Mil', c='darkred')

axs2[0].plot(P_out, output_array1[:, 1], label='v_out - Tal', c='darkblue', ls=':')
axs2[0].plot(P_out, output_array1[:, 2], label='sound_speed - Tal', c='darkred', ls=':')
axs2[0].legend()
axs2[0].grid()

axs2[1].plot(P_out, output_array[:, 0], label='m_dot - Mil', c='black')
axs2[1].plot(P_out, output_array1[:, 0], label='m_dot - Tal', c='black', ls=':')
axs2[1].legend()
axs2[1].grid()

plt.show()