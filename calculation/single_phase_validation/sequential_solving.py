# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.single_phase import SPStator, SPRotor, SPTeslaGeometry, SPTeslaOptions
from main_code.base_classes import BaseTeslaTurbine, CALCULATION_FOLDER
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os.path


# %%------------       INPUT DATA                         -----------------------------------------------------------> #
curr_geometry = SPTeslaGeometry()
curr_options = SPTeslaOptions()
curr_options.stator.iterate_phi = True
curr_geometry.stator.d_int = 0.2        # [m]
curr_options.rotor.integr_variable = 0.03
curr_options.rotor.n_rotor = 400

P_in = 22000000  # [Pa]
T_in = 423.15      # K
P_out = 21000000   # [Pa]

tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStator, rotor=SPRotor)
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)

tt.static_points[1].set_variable("P", P_out)

tt.geometry.d_main = 0.2
tt.geometry.throat_width = 0.0004934
tt.geometry.H_s = 0.00094
tt.geometry.alpha1 = 85

tt.rotor.rpm = 6000

tt.geometry.rotor.d_ratio = 2.5
tt.geometry.rotor.n_channels = 2


# %%------------             IMPORT EES RESULTS                ------------------------------------------------------> #
EES_FILE_FOLDER = os.path.join(CALCULATION_FOLDER, "single_phase_validation", "EES Results")
excel_file = os.path.join(EES_FILE_FOLDER, "rotor_array.xlsx")
df_ees = pd.read_excel(excel_file)
EES_array = df_ees.to_numpy()
EES_array = np.flip(EES_array, axis=0)

# %%------------     STATOR SOLVE                        -----------------------------------------------------------> #
tt.stator.solve()

# %%------------      ROTOR SOLVE                        -----------------------------------------------------------> #
tt.rotor.solve()
rotor_array = tt.rotor.get_rotor_array()
tt.evaluate_performances()
rpm = tt.rotor.rpm
P_out = rotor_array[-1, 11]
power = tt.power

# %%------------      EVALUATE ERROR                      -----------------------------------------------------------> #
abs_errors = np.zeros((len(rotor_array[:, 0]), 11))
rel_errors = np.zeros((len(rotor_array[:, 0]), 11))
r_old = rotor_array[0, 1]

h_0_ees = np.interp(r_old, EES_array[:, 0], EES_array[:, 2])
h_0_py = rotor_array[0, 12]

s_0_ees = np.interp(r_old, EES_array[:, 0], EES_array[:, 3])
s_0_py = rotor_array[0, 14]

p_ees_arr = np.interp(rotor_array[:, 1], EES_array[:, 0], EES_array[:, 1])

for i in range(len(rotor_array[:, 0])):

    if i > 0:

        r_curr = rotor_array[i, 1]

        # Pressures
        dp_py = rotor_array[i, 11] - rotor_array[i-1, 11]
        dp_ees = p_ees_arr[i] - p_ees_arr[i-1]

        abs_errors[i, 1] = (rotor_array[i, 11] - p_ees_arr[i])
        abs_errors[i, 2] = (dp_py - dp_ees)

        rel_errors[i, 1] = abs_errors[i, 1] / p_ees_arr[i]
        rel_errors[i, 2] = (dp_py - dp_ees) / dp_ees

        # Enthalpies
        h_ees = np.interp(r_curr, EES_array[:, 0], EES_array[:, 2]) - h_0_ees
        h_ees_old = np.interp(r_old, EES_array[:, 0], EES_array[:, 2]) - h_0_ees

        h_py = rotor_array[i, 12] - h_0_py
        dh_py = rotor_array[i, 12] - rotor_array[i - 1, 12]
        dh_ees = h_ees - h_ees_old

        abs_errors[i, 3] = (h_py - h_ees)
        abs_errors[i, 4] = (dh_py - dh_ees)

        rel_errors[i, 3] = (h_py - h_ees) / h_ees
        rel_errors[i, 4] = (dh_py - dh_ees) / dh_ees

        r_old = r_curr

abs_errors[:, 0] = rotor_array[:, 1]
rel_errors[:, 0] = rotor_array[:, 1]
abs_errors_df = pd.DataFrame(abs_errors)
rel_errors_df = pd.DataFrame(rel_errors)


# %%------------             TRAJECTORY PLOT                        -------------------------------------------------> #
fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

ax.plot(rotor_array[:, 24] *np.pi / 180, rotor_array[:, 1])
ax.set_rmax(0.1)
ax.set_rticks([0.025, 0.05, 0.075, 0.1])
ax.set_rlabel_position(-22.5)
ax.grid(True)

ax.set_title("theta_rel")

plt.show()

# %%------------             VELOCITY PLOT                        ---------------------------------------------------> #

plt.rcParams['figure.constrained_layout.use'] = True
fig1, (ax3, ax4) = plt.subplots(1, 2, figsize=(10, 6))

# Absolute Velocities
ax3.plot(rotor_array[:, 1], rotor_array[:, 3], label="Tangential Speed vt")
ax3.plot(rotor_array[:, 1], rotor_array[:, 4], label="Radial Speed vr")

ax3.set_xlabel("Radial Distance [m]", color = "Black", fontsize = 14)
ax3.set_ylabel("Velocity [m/s]", color = "Black", fontsize = 14)
ax3.set_xlim(0.04, 0.1)
ax3.set_ylim(0, 75)
ax3.legend(edgecolor = "Black", fontsize = 12)
ax3.grid()

# Relative Velocities
ax4.plot(rotor_array[:, 1], rotor_array[:, 5], label="Tangential Speed wt")
ax4.plot(rotor_array[:, 1], rotor_array[:, 6], label="Radial Speed wr")

ax4.set_xlabel("Radial Distance [m]", color = "Black", fontsize = 14)
ax4.set_ylabel("Velocity [m/s]", color = "Black", fontsize = 14)
ax4.set_xlim(0.04, 0.1)
ax4.set_ylim(0, 25)
ax4.legend(edgecolor = "Black", fontsize = 12)
ax4.grid()

plt.show()

# %%------------             VELOCITY PLOT 2                      ---------------------------------------------------> #

fig2 = plt.subplots()

plt.plot(rotor_array[:, 1], rotor_array[:, 2], label="Rotating Speed", c="Darkred")
plt.plot(rotor_array[:, 1], rotor_array[:, 9], label="Absolute Speed", c="Darkblue")
plt.plot(rotor_array[:, 1], rotor_array[:, 10], label="Relative Speed", c="Darkgreen")

plt.xlabel("Radial Distance [m]", color = "Black", fontsize = 14)
plt.ylabel("Velocity [m/s]", color = "Black", fontsize = 14)
plt.xlim(0.04, 0.1)
plt.ylim(0, 75)
plt.legend(edgecolor="Black", fontsize = 12)
plt.grid()

plt.show()

# %%------------             THERMODYNAMIC PLOT                   ---------------------------------------------------> #

fig3 = plt.subplots()

plt.plot(rotor_array[:, 1], rotor_array[:, 11]/rotor_array[0, 11], label="Pressure", c="Darkred")
plt.plot(rotor_array[3:, 1], rotor_array[3:, 13]/rotor_array[3, 13], label="Temperature", c="Darkgreen")
plt.plot(rotor_array[3:, 1], rotor_array[3:, 15]/rotor_array[3, 15], label="Density", c="Orange")
plt.plot(rotor_array[3:, 1], rotor_array[3:, 12]/rotor_array[3, 12], label="Enthalpy", c="Darkblue")
plt.plot(rotor_array[3:, 1], rotor_array[3:, 14]/rotor_array[3, 14], label="Entropy", c="Olive")


plt.ylabel("Non-Dimensional Parameter [-]", fontsize = 14)
plt.xlabel("Radial Distance [m]", fontsize = 14)
plt.xlim(0.04, 0.1)
plt.ylim(0.93, 1.02)
plt.legend(edgecolor = "Black", fontsize = 12)
plt.grid()

plt.show()


# %%------------            SAVE RESULTS                       ------------------------------------------------------> #
py_df = pd.DataFrame(rotor_array)
RESULT_FOLDER = os.path.join(CALCULATION_FOLDER, "single_phase_validation", "results")
excel_file = os.path.join(RESULT_FOLDER, "rotor_array.xlsx")

with pd.ExcelWriter(excel_file) as writer:

    py_df.to_excel(writer, sheet_name="Data")
    abs_errors_df.to_excel(writer, sheet_name="abs error")
    rel_errors_df.to_excel(writer, sheet_name="rel error")
