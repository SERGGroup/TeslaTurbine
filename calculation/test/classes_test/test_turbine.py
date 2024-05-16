# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.single_phase import SPStator, SPRotor, SPTeslaGeometry, SPTeslaOptions
from main_code.base_classes import BaseTeslaTurbine
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
import numpy as np


# %%------------       INPUT DATA                         -----------------------------------------------------------> #
curr_geometry = SPTeslaGeometry()
curr_options = SPTeslaOptions()
curr_options.stator.iterate_phi = True
curr_geometry.stator.d_int = 0.2        # [m]

tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStator, rotor=SPRotor)

P_in = 22000000  # [Pa]
P_out = np.linspace(20000000, 16000000, 20)
T_in = 423.15      # K
tt.T_in = T_in
tt.P_in = P_in
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)


tt.geometry.d_main = 0.2
tt.geometry.H_s = 0.00094
tt.geometry.alpha_stat = 85

# Design Parameters
tt.geometry.rotor.d_ratio = 2.5
tt.geometry.rotor.b_channel = 0.00007
tt.geometry.throat_width = 0.0004934
tt.rotor.dv_perc = 0.1315
tt.geometry.rotor.n_channels = 2

output_array = np.empty((len(P_out), 6))


# %%------------             CALCULATIONS                -----------------------------------------------------------> #
for i in tqdm(range(len(P_out))):

    tt.P_out = P_out[i]

    tt.iterate_pressure()
    rotor_array = tt.rotor.get_rotor_array()
    tt.evaluate_performances()

    output_array[i, 0] = P_out[i]
    output_array[i, 1] = tt.eta_tt
    output_array[i, 2] = tt.work
    output_array[i, 3] = tt.power
    output_array[i, 4] = tt.rotor.rpm
    output_array[i, 5] = tt.stator.m_dot_s


# %%------------             PLOT RESULTS                -----------------------------------------------------------> #
df = pd.DataFrame(output_array)
df.to_excel(excel_writer="test.xlsx")
