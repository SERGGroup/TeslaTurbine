# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
import matplotlib.pyplot as plt
from main_code.sub_classes.multi_phase import TPStator0D, TPRotor, TPTeslaGeometry, TPTeslaOptions, TPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
import pandas as pd

# %%------------   SETUP DATA                             -----------------------------------------------------------> #
curr_geometry = TPTeslaGeometry()
curr_geometry.stator.d_int = 0.5    # [m]
curr_geometry.stator.H_s = 0.0005   # [m]
curr_geometry.H_s = 0.0005          # [m]              TODO This needs to be cleaned up
curr_geometry.throat_width = 0.003  # [m]

curr_options = TPTeslaOptions()
curr_options.rotor.profile_rotor = True
curr_options.rotor.tp_epsilon_model = "chisholm"

P_in = 997233  # [Pa]
P_out = 652161

INLET_T = np.linspace(315.15, 323.15, 50)

tt1 = BaseTeslaTurbine("R1234ze", curr_geometry, curr_options, stator=TPStatorMil, rotor=TPRotor)

tt1.points[0].set_variable("P", P_in)
tt1.static_points[1].set_variable("P", P_out)

tt1.P_in = P_in                                     # TODO This needs to be cleaned up too
tt1.stator.stator_eff = 0.81

j = 7
output_array1 = np.empty((len(INLET_T), j))

# %%------------   COMPARE STATOR MODELS                     -------------------------------------------------------> #
for i in range(len(INLET_T)):

    tt1.points[0].set_variable("T", INLET_T[i])
    tt1.T_in = INLET_T[i]                           # TODO This needs to be cleaned up too

    tt1.stator.solve()

    output_array1[i, 0] = INLET_T[i]

    output_array1[i, 1] = tt1.stator.m_dot_s
    output_array1[i, 2] = tt1.stator.stator_eff
    output_array1[i, 3] = tt1.stator.out_speed
    output_array1[i, 4] = tt1.stator.x_out
    output_array1[i, 5] = tt1.stator.x_in
    output_array1[i, 6] = tt1.stator.phi_n

df = pd.DataFrame(output_array1, columns=["inlet T", "m_dot_s", "stator_eff", "out_speed", "x_out", "x_in", "phi_n"])

# %%------------   PLOT COMPARISON                           -------------------------------------------------------> #
fig = plt.subplots()

plt.plot(output_array1[:, 0], output_array1[:, 4], c='darkblue')
plt.grid()

plt.show()
