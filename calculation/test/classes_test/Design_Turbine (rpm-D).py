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
tt.geometry.alpha_stat = 85

# Design Parameters
tt.geometry.rotor.d_ratio = 2.5
tt.geometry.rotor.b_channel = 0.00007
tt.geometry.throat_width = 0.0005

D_main = np.linspace(0.2, 0.5, 10)
dw_perc = np.linspace(-0.3, 0.3, 20)

output_array = np.empty((len(D_main) * len(dw_perc), 6))

# %%------------             CALCULATIONS                -----------------------------------------------------------> #
for i in tqdm(range(len(D_main))):

    tt.geometry.d_main = D_main[i]
    curr_geometry.stator.d_int = D_main[i]

    for j in range(len(dw_perc)):

        tt.rotor.dv_perc = dw_perc[j]

        tt.iterate_pressure()
        rotor_array = tt.rotor.get_rotor_array()
        tt.evaluate_performances()

        output_array[len(dw_perc) * i + j, 0] = D_main[i]
        output_array[len(dw_perc) * i + j, 1] = dw_perc[j]
        output_array[len(dw_perc) * i + j, 2] = tt.Eta_tesla_ss
        output_array[len(dw_perc) * i + j, 3] = tt.work
        output_array[len(dw_perc) * i + j, 4] = tt.power
        output_array[len(dw_perc) * i + j, 5] = tt.rotor.rpm

# %%------------             PLOT RESULTS                -----------------------------------------------------------> #
print(output_array)
