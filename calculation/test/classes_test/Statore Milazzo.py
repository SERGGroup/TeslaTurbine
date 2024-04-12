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
P_out = np.linspace(18000000, 4000000, 50)
T_in = 423.15      # K
tt.T_in = T_in
tt.P_in = P_in
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)

tt.geometry.d_main = 0.2
tt.geometry.H_s = 0.00094
tt.geometry.alpha1 = 85
tt.stator.stator_eff = 0.9409

# Design Parameters
tt.geometry.rotor.d_ratio = 2.5
tt.geometry.rotor.b_channel = 0.00007
tt.geometry.throat_width = 0.0004934
tt.rotor.rpm = 6000
tt.geometry.rotor.n_channels = 2

output_array = np.empty((len(P_out), 10))


# %%------------             CALCULATIONS                -----------------------------------------------------------> #
for i in tqdm(range(len(P_out))):

    tt.P_out = P_out[i]

    tt.iterate_pressure()
    rotor_array = tt.rotor.get_rotor_array()
    tt.evaluate_performances()

    output_array[i, 0] = P_out[i]
    output_array[i, 1] = tt.Eta_tesla_ss
    output_array[i, 2] = tt.work
    output_array[i, 3] = tt.power
    output_array[i, 4] = tt.rotor.rpm
    output_array[i, 5] = tt.stator.m_dot_s
    output_array[i, 6] = tt.stator.out_speed
    output_array[i, 7] = tt.points[1].get_variable("P")
    output_array[i, 8] = tt.points[1].get_variable("P") - tt.points[2].get_variable("P")
    output_array[i, 9] = tt.stator.Ma_1

py_df = pd.DataFrame(output_array, columns=['P_out', 'Eta_Tesla_ss', 'Work', 'Power', 'RPM', 'm_dot',
                                            'Stator_Outlet_Speed', 'P[1]', 'DP_Gap', 'Ma[1]'])

# %%------------       OLD STATOR                         -----------------------------------------------------------> #

curr_geometry1 = SPTeslaGeometry()
curr_options1 = SPTeslaOptions()
curr_options1.stator.iterate_phi = True
curr_geometry1.stator.d_int = 0.2        # [m]
curr_options1.rotor.integr_variable = 0.03


tt1 = BaseTeslaTurbine("CarbonDioxide", curr_geometry1, curr_options1, stator=SPStator, rotor=SPRotor)

P_in1 = 22000000  # [Pa]
P_out1 = np.linspace(18000000, 4000000, 50)
T_in1 = 423.15      # K
tt1.T_in = T_in1
tt1.P_in = P_in1
tt1.points[0].set_variable("P", P_in1)
tt1.points[0].set_variable("T", T_in1)

tt1.geometry.d_main = 0.2
tt1.geometry.H_s = 0.00094
tt1.geometry.alpha1 = 85
tt1.rotor.rpm = 6000

# Design Parameters
tt1.geometry.rotor.d_ratio = 2.5
tt1.geometry.rotor.b_channel = 0.00007
tt1.geometry.throat_width = 0.0004934
tt1.geometry.rotor.n_channels = 2

output_array1 = np.empty((len(P_out1), 10))


# %%------------             CALCULATIONS                -----------------------------------------------------------> #
for i in tqdm(range(len(P_out1))):

    tt1.P_out = P_out1[i]

    tt1.iterate_pressure()
    rotor_array1 = tt1.rotor.get_rotor_array()
    tt1.evaluate_performances()

    output_array1[i, 0] = P_out1[i]
    output_array1[i, 1] = tt1.Eta_tesla_ss
    output_array1[i, 2] = tt1.work
    output_array1[i, 3] = tt1.power
    output_array1[i, 4] = tt1.rotor.rpm
    output_array1[i, 5] = tt1.stator.m_dot_s
    output_array1[i, 6] = tt1.stator.v_out
    output_array1[i, 7] = tt1.points[1].get_variable("P")
    output_array1[i, 8] = tt1.points[1].get_variable("P") - tt1.points[2].get_variable("P")
    output_array1[i, 9] = tt1.stator.Ma_1

py_df1 = pd.DataFrame(output_array1, columns=['P_out', 'Eta_Tesla_ss', 'Work', 'Power', 'RPM', 'm_dot',
                                            'Stator_Outlet_Speed', 'P[1]', 'DP_Gap', 'Ma[1]'])