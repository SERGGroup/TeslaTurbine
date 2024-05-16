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
curr_geometry.stator.d_int = 0.2        # [m]
curr_options.rotor.integr_variable = 0.03


tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStator, rotor=SPRotor)

P_in = 22000000  # [Pa]
P_out = np.linspace(19139156, 14524628, 20)
T_in = 423.15      # K
tt.T_in = T_in
tt.P_in = P_in
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)

tt.geometry.d_main = 0.2
tt.geometry.H_s = 0.00094
tt.geometry.alpha_stat = 85
tt.rotor.gap_losses_control = True

# Design Parameters
tt.geometry.rotor.d_ratio = 2.5
tt.geometry.rotor.b_channel = 0.00007
tt.geometry.throat_width = 0.0004934
tt.geometry.rotor.n_channels = 2
tt.rotor.rpm = 6000

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

py_df = pd.DataFrame(output_array, columns=['P_out', 'Eta_Tesla_ss', 'Work', 'Power', 'RPM', 'm_dot'])

# %%------------             IMPORT EES RESULTS                ------------------------------------------------------> #

df1 = pd.read_excel("C:/Users/iggig/Desktop/ees_test.xlsx")

# Extracting Array from Dataframe
EES_array = df1.to_numpy()

# %%------------------             PLOT RESULTS                ------------------------------------------------------> #
plt.rcParams['figure.constrained_layout.use'] = True
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10, 6))


ax1.plot(output_array[:, 0], output_array[:, 1], label = 'Python Model', color = "Darkred")
ax1.plot(EES_array[:, 7], EES_array[:, 3], label = 'EES Model', color = "Darkblue")

ax1.set_xlabel("Outlet Pressure [Pa]", color="Black", fontsize = 14)
ax1.set_ylabel("Eta_ss [-]", color="Black", fontsize = 14)
ax1.set_xlim(14200000, 19500000)
ax1.set_ylim(0.38, 0.55)
ax1.legend(edgecolor = "Black", fontsize = 12)
ax1.grid()

ax2.plot(output_array[:, 0], output_array[:, 3], label = 'Python Model', color = "Darkred")
ax2.plot(EES_array[:, 7], EES_array[:, 1], label = 'EES Model', color = "Darkblue")

ax2.set_xlabel("Outlet Pressure [Pa]", color = "Black", fontsize = 14)
ax2.set_ylabel("Power [W]", color = "Black", fontsize = 14)
ax2.set_xlim(14200000, 19500000)
ax2.set_ylim(75, 240)
ax2.legend(edgecolor = "Black", fontsize = 12)
ax2.grid()

plt.show()
