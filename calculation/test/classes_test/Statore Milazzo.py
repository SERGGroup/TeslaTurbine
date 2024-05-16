# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
import numpy as np
from main_code.sub_classes.single_phase import SPStator, SPRotor, SPTeslaGeometry, SPTeslaOptions, SPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd


# %%------------       MILAZZO STATOR                         -------------------------------------------------------> #

curr_geometry = SPTeslaGeometry()
curr_options = SPTeslaOptions()
curr_options.stator.iterate_phi = True
curr_geometry.stator.d_int = 0.2        # [m]
curr_options.rotor.integr_variable = 0.03

tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStatorMil, rotor=SPRotor)

P_in = 22000000  # [Pa]
P_out = np.linspace(21000000, 10000000, 200)
T_in = 423.15      # K
tt.T_in = T_in
tt.P_in = P_in
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)
tt.stator.stator_eff = 0.92

tt.geometry.d_main = 0.2
tt.geometry.H_s = 0.00094
tt.geometry.alpha_stat = 85

# Design Parameters
tt.geometry.rotor.d_ratio = 2.5
tt.geometry.rotor.b_channel = 0.00007
tt.geometry.throat_width = 0.0004934
tt.rotor.rpm = 6000
tt.geometry.rotor.n_channels = 2

# %%------------       OLD STATOR                         -----------------------------------------------------------> #
curr_geometry1 = SPTeslaGeometry()
curr_options1 = SPTeslaOptions()
curr_options1.stator.iterate_phi = True
curr_geometry1.stator.d_int = 0.2        # [m]
curr_options1.rotor.integr_variable = 0.03

tt1 = BaseTeslaTurbine("CarbonDioxide", curr_geometry1, curr_options1, stator=SPStator, rotor=SPRotor)

tt1.T_in = T_in
tt1.P_in = P_in
tt1.points[0].set_variable("P", P_in)
tt1.points[0].set_variable("T", T_in)

tt1.geometry.d_main = tt.geometry.d_main
tt1.geometry.H_s = tt.geometry.H_s
tt1.geometry.alpha_stat = tt.geometry.alpha_stat
tt1.rotor.rpm = tt.rotor.rpm

# Design Parameters
tt1.geometry.rotor.d_ratio = tt.geometry.rotor.d_ratio
tt1.geometry.rotor.b_channel = tt.geometry.rotor.b_channel
tt1.geometry.throat_width = tt.geometry.throat_width
tt1.geometry.rotor.n_channels = tt.geometry.rotor.n_channels


# %%------------             CALCULATIONS                -----------------------------------------------------------> #
output_array1 = np.empty((len(P_out), 10))

for i in tqdm(range(len(P_out))):

    tt1.P_out = P_out[i]

    tt1.solve_with_stator_outlet_pressure(P_out[i])
    rotor_array1 = tt1.rotor.get_rotor_array()
    tt1.evaluate_performances()

    output_array1[i, 0] = P_out[i]
    output_array1[i, 1] = tt1.eta_tt
    output_array1[i, 2] = tt1.work
    output_array1[i, 3] = tt1.power
    output_array1[i, 4] = tt1.rotor.rpm
    output_array1[i, 5] = tt1.stator.m_dot_s
    output_array1[i, 6] = tt1.stator.v_out
    output_array1[i, 7] = tt1.points[1].get_variable("P")
    output_array1[i, 8] = tt1.points[1].get_variable("P") - tt1.points[2].get_variable("P")
    output_array1[i, 9] = tt1.static_points[1].get_variable("P")

py_df1 = pd.DataFrame(output_array1, columns=['P_out', 'Eta_Tesla_ss', 'Work', 'Power', 'RPM', 'm_dot',
                                            'Stator_Outlet_Speed', 'P[1]', 'DP_Gap', 'Ma[1]'])

# %%------------             CALCULATIONS                -----------------------------------------------------------> #
output_array = np.empty((len(P_out), 10))

for i in tqdm(range(len(P_out))):

    tt.stator.p_out = P_out[i]
    # tt.stator.stator_eff = output_array1[i, 9]

    tt.solve_with_stator_outlet_pressure(P_out[i])
    rotor_array = tt.rotor.get_rotor_array()
    tt.evaluate_performances()

    output_array[i, 0] = P_out[i]
    output_array[i, 1] = tt.eta_tt
    output_array[i, 2] = tt.work
    output_array[i, 3] = tt.power
    output_array[i, 4] = tt.rotor.rpm
    output_array[i, 5] = tt.stator.m_dot_s
    output_array[i, 6] = tt.stator.speed_out.v
    output_array[i, 7] = tt.points[1].get_variable("P")
    output_array[i, 8] = tt.points[1].get_variable("P") - tt.points[2].get_variable("P")
    output_array[i, 9] = tt.static_points[1].get_variable("P")

py_df = pd.DataFrame(output_array, columns=['P_out', 'Eta_Tesla_ss', 'Work', 'Power', 'RPM', 'm_dot',
                                            'Stator_Outlet_Speed', 'P[1]', 'DP_Gap', 'Ma[1]'])


# %%------------             PLOT                        -----------------------------------------------------------> #
i = 0
j = 5
fig = plt.subplots()

plt.plot(output_array[:, i], output_array[:, j], color='blue')
plt.plot(output_array1[:, i], output_array1[:, j], color='red')

plt.legend(['MILAZZO', 'OLD'], loc='lower left')
plt.show()
