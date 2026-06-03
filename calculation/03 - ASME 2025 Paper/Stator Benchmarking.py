# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
import numpy as np
from main_code.sub_classes.single_phase import SPStator, SPRotor, SPTeslaGeometry, SPTeslaOptions, SPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd

# %%------------       INPUT DATA                         -----------------------------------------------------------> #
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


P_in = 10100000  # [Pa]
T_in_c = 94  # [°C]
T_in = T_in_c + 273.15  # [K]

m_rif = 0.157

P_out = np.linspace(10000000, 4000000, 100)

tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStator, rotor=SPRotor)
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)

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

tt.stator.stator_eff = 0.9
tt.rotor.gap_losses_control = True
tt.include_extra_losses = False

tt1 = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStatorMil, rotor=SPRotor)
tt1.points[0].set_variable("P", P_in)
tt1.points[0].set_variable("T", T_in)

tt1.rotor.rpm = 15000

tt1.geometry.stator.Z_stat = Z_stat

tt1.geometry.throat_width = throat_width

tt1.geometry.disc_thickness = 0.0001
tt1.geometry.rotor.n_discs = 1
tt1.geometry.H_s = 0.00015
tt1.m_dot_tot = m_rif
tt1.geometry.n_channels = tt.n_packs

tt1.geometry.alpha1 = alpha_in

tt1.stator.stator_eff = 0.95
tt1.rotor.gap_losses_control = True
tt1.include_extra_losses = False

output_array = np.empty((len(P_out), 5))

tt1.P_in = P_in
tt1.P_out = P_out
tt1.T_in = T_in

# %%------------             CALCULATIONS                -----------------------------------------------------------> #
for i in tqdm(range(len(P_out))):

    tt.P_out = P_out[i]
    tt1.P_out = P_out[i]

    tt.static_points[1].set_variable("P", P_out[i])
    tt1.static_points[1].set_variable("P", P_out[i])

    tt.stator.solve()
    tt1.stator.solve()

    output_array[i, 0] = P_out[i]
    output_array[i, 1] = tt.stator.eta_stat
    output_array[i, 2] = tt.stator.m_dot_s
    output_array[i, 3] = tt1.stator.eta_stat
    output_array[i, 4] = tt1.stator.m_dot_s


# %%------------             COMPARISON WITH EES                        --------------------------------------------> #
P = np.linspace(19500000, 21500000, 21)
Eta = [0.9347, 0.9344, 0.9341, 0.9339, 0.9336, 0.9332, 0.933, 0.9327, 0.9322, 0.9319, 0.9314, 0.931, 0.9306, 0.9301, 0.9294, 0.9287, 0.9281, 0.9273, 0.9264, 0.9253, 0.924]
m_dot = [0.07242, 0.07112, 0.06977, 0.06839, 0.06696, 0.06548, 0.06396, 0.06239, 0.06075, 0.05906, 0.0573, 0.05546, 0.05355, 0.05155, 0.04944, 0.04723, 0.04488, 0.04239, 0.03971, 0.03682, 0.03366]

# %%------------             PLOT                        -----------------------------------------------------------> #

plt.plot(output_array[:, 0] / 100000, output_array[:, 2], color="darkblue")
plt.plot(output_array[:, 0] / 100000, output_array[:, 4], color="darkred")
plt.xlabel("Stator Outlet Pressure [bar]", fontsize=12)
plt.ylabel("Mass Flow Rate [kg/s]", fontsize=12)
plt.grid()

plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.tight_layout()
plt.legend(["[29] Stator Model", "Current Work"], loc="lower left", fontsize=11, edgecolor="black")

plt.show()
