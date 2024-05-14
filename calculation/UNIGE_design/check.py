# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
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
curr_geometry.stator.d_int = 0.3        # [m]
curr_options.rotor.integr_variable = 0.03


tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStator, rotor=SPRotor)
tt.rotor.gap_losses_control = False

P_in = 10000000     # [Pa]
P_out = 3800000     # [Pa]
T_in = 367.15       # [K]

tt.points[1].set_variable("P", P_in)
tt.points[1].set_variable("T", T_in)

tt.static_points[1].set_variable("P", P_in)
tt.static_points[1].set_variable("T", T_in)

ss = tt.points[1].get_variable("c")
alpha = 87.
alpha_rad = alpha * np.pi / 180
M = 0.8
tt.geometry.rotor.n_channels = 100
tt.stator.m_dot_s = 0.1

v_in = M * ss
tt.stator.speed_out.init_from_codes("v", v_in, "alpha", alpha_rad)

tt.geometry.d_main = 0.3

tt.rotor.rpm = 5000

# Design Parameters
D_ratio = np.linspace(2.5, 4.5, 10)
tt.geometry.rotor.b_channel = 0.0001

output_array = np.empty((len(D_ratio), 6))


# %%------------             CALCULATIONS                -----------------------------------------------------------> #
for i in tqdm(range(len(D_ratio))):

    tt.geometry.rotor.d_ratio = D_ratio[i]

    tt.rotor.solve()
    rotor_array = tt.rotor.get_rotor_array()
    tt.evaluate_rotor_performances()

    output_array[i, 0] = D_ratio[i]
    output_array[i, 1] = tt.rotor.eta_tt
    output_array[i, 2] = tt.work
    output_array[i, 3] = tt.power
    output_array[i, 4] = tt.rotor.rpm
    output_array[i, 5] = tt.rotor.m_dot_ch

py_df = pd.DataFrame(output_array, columns=['D_ratio', 'Eta_Tesla_ss', 'Work', 'Power', 'RPM', 'm_dot'])

# %%------------------             PLOT RESULTS                ------------------------------------------------------> #
