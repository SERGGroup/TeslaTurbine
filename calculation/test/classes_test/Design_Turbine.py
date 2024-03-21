# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
import numpy as np
from main_code.sub_classes.single_phase import SPStator, SPRotor, SPTeslaGeometry, SPTeslaOptions
from main_code.base_classes import BaseTeslaTurbine
import matplotlib.pyplot as plt
from tqdm import tqdm

# %%------------       INPUT DATA                         -----------------------------------------------------------> #
curr_geometry = SPTeslaGeometry()
curr_options = SPTeslaOptions()
curr_options.stator.iterate_phi = True
curr_geometry.stator.d_int = 0.2        # [m]

tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStator, rotor=SPRotor)

P_in = 10000000                               # [Pa]
tt.P_in = P_in
T_in = 367.05                                 # [K]
tt.T_in = T_in
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)
tt.P_out = 3800000                            # [Pa]

tt.geometry.rotor.n_channels = 2
tt.geometry.H_s = 0.001
tt.geometry.alpha1 = 85

# Design Parameters
tt.geometry.rotor.d_ratio = 2.5
tt.geometry.d_main = 0.2
tt.geometry.rotor.b_channel = 0.00007
tt.geometry.throat_width = 0.0004934
tt.rotor.dv_perc = 0.1315

tt.iterate_pressure()
rotor_array = tt.rotor.get_rotor_array()

# output_array = np.empty((len(P_in), 4))

# %%------------             CALCULATIONS                -----------------------------------------------------------> #
# for i in tqdm(range(len(P_in))):

#    tt.P_in = P_in[i]
#    tt.points[0].set_variable("P", P_in[i])
#    tt.points[0].set_variable("T", T_in)

#   tt.iterate_pressure()
#   rotor_array = tt.rotor.get_rotor_array()

#    output_array[i, 0] = tt.stator.p_out
#    output_array[i, 1] = tt.stator.m_dot_s
#    output_array[i, 2] = rotor_array[-1, 9]
#    output_array[i, 3] = rotor_array[-1, 13]

# %%------------             PLOT RESULTS                -----------------------------------------------------------> #
fig, ax = plt.subplots()

# plt.plot(P_in, output_array[:, 1])
plt.show()
# print(output_array)
