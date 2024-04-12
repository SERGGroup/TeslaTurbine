# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
import numpy as np
from main_code.sub_classes.single_phase import SPStator, SPRotor, SPTeslaGeometry, SPTeslaOptions
from main_code.base_classes import BaseTeslaTurbine
import matplotlib.pyplot as plt
import tqdm as tqdm
import pandas as pd

# %%------------       INPUT DATA                         -----------------------------------------------------------> #
curr_geometry = SPTeslaGeometry()
curr_options = SPTeslaOptions()
curr_options.stator.iterate_phi = True
curr_geometry.stator.d_int = 0.2        # [m]

P_in = 22000000  # [Pa]
T_in = 423.15      # K
# P_out = 21000000   # [Pa]
P_out = np.linspace(19500000, 21500000, 21)

tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStator, rotor=SPRotor)
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)

tt.geometry.d_main = 0.2
tt.geometry.throat_width = 0.0004934
tt.geometry.H_s = 0.00094
tt.geometry.alpha1 = 85

tt.geometry.rotor.d_ratio = 2.5
tt.geometry.rotor.n_channels = 2

output_array = np.empty((len(P_out), 3))

# %%------------             CALCULATIONS                -----------------------------------------------------------> #
for i in range(len(P_out)):

    tt.P_out = P_out[i]

    tt.static_points[1].set_variable("P", P_out[i])
    tt.stator.solve()

    output_array[i, 0] = P_out[i]
    output_array[i, 1] = tt.stator.eta_stat
    output_array[i, 2] = tt.stator.m_dot_s

# %%------------             COMPARISON WITH EES                        --------------------------------------------> #

P = np.linspace(19500000, 21500000, 21)
Eta = [0.9347, 0.9344, 0.9341, 0.9339, 0.9336, 0.9332, 0.933, 0.9327, 0.9322, 0.9319, 0.9314, 0.931, 0.9306, 0.9301, 0.9294, 0.9287, 0.9281, 0.9273, 0.9264, 0.9253, 0.924]
m_dot = [0.07242, 0.07112, 0.06977, 0.06839, 0.06696, 0.06548, 0.06396, 0.06239, 0.06075, 0.05906, 0.0573, 0.05546, 0.05355, 0.05155, 0.04944, 0.04723, 0.04488, 0.04239, 0.03971, 0.03682, 0.03366]

# %%------------             PLOT                        -----------------------------------------------------------> #
fig = plt.subplots()

plt.plot(output_array[:, 0], output_array[:,2], color = "black")
plt.plot(P, m_dot, color = "darkred")
plt.xlabel("Stator Outlet Pressure [Pa]", fontsize = 13)
plt.ylabel("Mass Flow Rate [kg/s]", fontsize = 13)
plt.xlim(19500000, 21500000)
plt.ylim(0.03, 0.075)
plt.grid()

plt.legend(["Python Model", "EES Model"], loc = "lower left", fontsize = 13, edgecolor = "black")

plt.show()
