# %%------------   IMPORT CLASSES
from main_code.sub_classes.single_phase import SPStator, SPRotor, SPTeslaGeometry, SPTeslaOptions
from main_code.base_classes import BaseTeslaTurbine, CALCULATION_FOLDER
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
import numpy as np
import os.path


# %%------------       INPUT DATA                         -----------------------------------------------------------> #
curr_geometry = SPTeslaGeometry()
curr_options = SPTeslaOptions()
curr_options.stator.iterate_phi = True
curr_geometry.stator.d_int = 0.2        # [m]
curr_options.rotor.integr_variable = 0.03

tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStator, rotor=SPRotor)

P_in = 22000000   # [Pa]
P_out = 21000000
T_in = 423.15  # K
tt.T_in = T_in
tt.P_in = P_in
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)

d_main = 0.2
tt.geometry.H_s = 0.00094
tt.geometry.alpha_stat = 85
rpm = np.linspace(4000, 7000, 10)

# Design Parameters
tt.geometry.rotor.d_ratio = 2.5
tt.geometry.rotor.b_channel = 0.00007
tt.geometry.throat_width = 0.0004934
tt.rotor.dv_perc = 0.1315
tt.geometry.rotor.n_channels = 2

output_array = np.empty((len(rpm), 9))


# %%------------             CALCULATIONS                -----------------------------------------------------------> #
for i in tqdm(range(len(rpm))):

    tt.rotor.rpm = rpm[i]

    tt.solve_with_stator_outlet_pressure(P_out)
    rotor_array = tt.rotor.get_rotor_array()
    tt.evaluate_performances()

    output_array[i, 0] = rpm[i]
    output_array[i, 1] = tt.Eta_tesla_ss
    output_array[i, 2] = tt.Eta_tesla_ss_dh
    output_array[i, 3] = tt.work
    output_array[i, 4] = tt.power
    output_array[i, 5] = tt.rotor.rpm
    output_array[i, 6] = tt.stator.m_dot_s
    output_array[i, 7] = tt.points[-1].get_variable("P")
    output_array[i, 8] = tt.points[-1].get_variable("T")

# %%------------   PLOT RESULTS                         -------------------------------------------------------------> #
fig2, ax2 = plt.subplots()
ax2.plot(output_array[:, 1], output_array[:, 2], c='darkblue', label='No-Meta')
plt.legend()
plt.grid()

plt.show()

# %%------------             IMPORT EES RESULTS                ------------------------------------------------------> #
py_df = pd.DataFrame(output_array, columns=['rpm', 'Eta_Tesla_ss', 'Work', 'Power', 'RPM', 'm_dot', "P_out", "T_out"])
RESULT_FOLDER = os.path.join(CALCULATION_FOLDER, "single_phase_validation", "results")
excel_file = os.path.join(RESULT_FOLDER, "rpm_variation_1.xlsx")
py_df.to_excel(excel_file)
