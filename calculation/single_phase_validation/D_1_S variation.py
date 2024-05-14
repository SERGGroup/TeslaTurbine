# %%------------   IMPORT CLASSES
import numpy as np
from main_code.sub_classes.single_phase import SPStator, SPRotor, SPTeslaGeometry, SPTeslaOptions
from main_code.base_classes import BaseTeslaTurbine
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd

# %%------------       INPUT DATA                         -----------------------------------------------------------> #
curr_geometry = SPTeslaGeometry()
curr_options = SPTeslaOptions()
curr_options.stator.iterate_phi = True
curr_geometry.stator.d_int = 0.2        # [m]
curr_options.rotor.integr_variable = 0.03


tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStator, rotor=SPRotor)

P_in = 22000000   # [Pa]
P_out = 20000000
T_in = 423.15  # K
tt.T_in = T_in
tt.P_in = P_in
tt.P_out = P_out
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)

d_main = np.linspace(0.1, 0.4, 50)
tt.geometry.H_s = 0.00094
tt.geometry.alpha1 = 85
tt.rotor.omega = 628.3185

# Design Parameters
tt.geometry.rotor.d_ratio = 2.5
tt.geometry.rotor.b_channel = 0.00007
tt.geometry.throat_width = 0.0004934
tt.rotor.dv_perc = 0.1315
tt.geometry.rotor.n_channels = 2

output_array = np.empty((len(d_main), 6))

# %%------------             CALCULATIONS                -----------------------------------------------------------> #
for i in tqdm(range(len(d_main))):

    tt.geometry.d_main = d_main[i]

    tt.iterate_pressure()
    rotor_array = tt.rotor.get_rotor_array()
    tt.evaluate_performances()

    output_array[i, 0] = d_main[i]
    output_array[i, 1] = tt.eta_tt
    output_array[i, 2] = tt.work
    output_array[i, 3] = tt.power
    output_array[i, 4] = tt.rotor.rpm
    output_array[i, 5] = tt.stator.m_dot_s

py_df = pd.DataFrame(output_array, columns=['d_main', 'Eta_Tesla_ss', 'Work', 'Power', 'RPM', 'm_dot'])

# %%------------             IMPORT EES RESULTS                ------------------------------------------------------> #

df1 = pd.read_excel("test.xlsx")

# Extracting Array from Dataframe
EES_array = df1.to_numpy()
