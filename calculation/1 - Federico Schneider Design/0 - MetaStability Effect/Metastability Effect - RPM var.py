# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.multi_phase import TPRotor, TPTeslaGeometry, TPTeslaOptions, TPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt


# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #
curr_geometry = TPTeslaGeometry()
curr_geometry.stator.d_int = 0.5            # [m]
curr_geometry.throat_width = 0.003          # [m]
curr_geometry.rotor.n_channels = 1          # [-]
curr_geometry.rotor.b_channel = 0.0003      # [m]
curr_geometry.rotor.roughness = 0.0000005   # [m]

curr_options1 = TPTeslaOptions()
curr_options1.rotor.profile_rotor = True
curr_options1.rotor.sp_check = False
curr_options1.rotor.tp_epsilon_model = "sarti"
curr_options1.stator.metastability_check = False

curr_options2 = TPTeslaOptions()
curr_options2.rotor.profile_rotor = True
curr_options2.rotor.sp_check = False
curr_options2.rotor.tp_epsilon_model = "sarti"
curr_options2.stator.metastability_check = True

fluid = "R1234ze(Z)"
tesla_turbine1 = BaseTeslaTurbine("R1234ze(Z)", curr_geometry, curr_options1, stator=TPStatorMil, rotor=TPRotor)
tesla_turbine2 = BaseTeslaTurbine("R1234ze(Z)", curr_geometry, curr_options2, stator=TPStatorMil, rotor=TPRotor)

P_in = 997233         # [Pa]
x_in = 0              # [-]
P_out = 427308        # [Pa]

tesla_turbine1.points[0].set_variable("P", P_in)
tesla_turbine1.points[0].set_variable("x", x_in)
tesla_turbine1.stator.stator_eff = 0.81
tesla_turbine1.rotor.gap_losses_control = True

tesla_turbine2.points[0].set_variable("P", P_in)
tesla_turbine2.points[0].set_variable("x", x_in)
tesla_turbine2.stator.stator_eff = 0.81
tesla_turbine2.rotor.gap_losses_control = True

dv_perc = np.linspace(-0.8, 0.3, 10)

output_array1 = np.empty((len(dv_perc), 9))
output_array2 = np.empty((len(dv_perc), 9))

Eta_meta = np.empty((len(dv_perc)))
Eta_meta_I = np.empty((len(dv_perc)))

Eta_nometa = np.empty((len(dv_perc)))
Eta_nometa_I = np.empty((len(dv_perc)))


# %%------------   TURBINE-NO_META SOLVING                         --------------------------------------------------> #
for i in tqdm(range(len(dv_perc))):

    tesla_turbine1.rotor.dv_perc = dv_perc[i]

    tesla_turbine1.P_in = P_in
    tesla_turbine1.P_out = P_out
    tesla_turbine1.iterate_pressure()
    tesla_turbine1.evaluate_performances()

    Eta_nometa[i] = tesla_turbine1.eta_tt
    Eta_nometa_I[i] = tesla_turbine1.eta_ss

    if Eta_nometa[i] > Eta_nometa[i-1] or i == 0:

        rotor_array = tesla_turbine1.rotor.get_rotor_array()

    output_array1[i, 0] = dv_perc[i]
    output_array1[i, 1] = tesla_turbine1.eta_tt
    output_array1[i, 2] = tesla_turbine1.work
    output_array1[i, 3] = tesla_turbine1.power
    output_array1[i, 4] = tesla_turbine1.rotor.rpm
    output_array1[i, 5] = tesla_turbine1.stator.m_dot_s
    output_array1[i, 6] = tesla_turbine1.points[1].get_variable("rho")
    output_array1[i, 7] = tesla_turbine1.static_points[1].get_variable("p")
    output_array1[i, 8] = tesla_turbine1.stator.speed_out.v

py_df = pd.DataFrame(output_array1, columns=['dv_perc', 'Eta_Tesla_ss', 'Work', 'Power', 'RPM', 'm_dot', 'rho_1', 'P_1', 'v'])


# %%------------   TURBINE-META SOLVING                         -----------------------------------------------------> #
for i in tqdm(range(len(dv_perc))):

    tesla_turbine2.rotor.dv_perc = dv_perc[i]

    tesla_turbine2.P_in = P_in
    tesla_turbine2.P_out = P_out
    tesla_turbine2.iterate_pressure()
    tesla_turbine2.evaluate_performances()

    Eta_meta[i] = tesla_turbine2.eta_tt
    Eta_meta_I[i] = tesla_turbine2.eta_ss

    if Eta_meta[i] > Eta_meta[i - 1] or i == 0:

        rotor_array2 = tesla_turbine2.rotor.get_rotor_array()

    output_array2[i, 0] = dv_perc[i]
    output_array2[i, 1] = tesla_turbine2.eta_tt
    output_array2[i, 2] = tesla_turbine2.work
    output_array2[i, 3] = tesla_turbine2.power
    output_array2[i, 4] = tesla_turbine2.rotor.rpm
    output_array2[i, 5] = tesla_turbine2.stator.m_dot_s
    output_array2[i, 6] = tesla_turbine2.points[1].get_variable("rho")
    output_array2[i, 7] = tesla_turbine2.static_points[1].get_variable("p")
    output_array2[i, 8] = tesla_turbine2.stator.speed_out.v


py_df2 = pd.DataFrame(output_array2, columns=['dv_perc', 'Eta_Tesla_ss', 'Work', 'Power', 'RPM', 'm_dot', 'rho_1', 'P_1', 'v'])


# %%------------   PLOT RESULTS                         -------------------------------------------------------------> #
param = 1
fig, ax = plt.subplots()
ax.plot(output_array2[:, 0], output_array2[:, param], c='darkblue', label='Meta')
ax.plot(output_array1[:, 0], output_array1[:, param], c='darkred', label='No-Meta')
ax.set(xlabel='dv_perc [-]', ylabel='Velocity [m/s]')

plt.legend()
plt.grid()

plt.show()


# %%------------   PLOT RESULTS                         -------------------------------------------------------------> #
fig1, axs = plt.subplots(1, 3, constrained_layout=True)

axs[0].plot(output_array1[:, 0], output_array1[:, 6], c='darkblue', label='No-Meta')
axs[0].plot(output_array2[:, 0], output_array2[:, 6], c='darkred', label='Meta')
axs[0].set(xlabel='dv_perc [-]', ylabel='Outlet Density [kg/m3]')
axs[0].legend(loc='upper right')

axs[1].plot(output_array1[:, 0], output_array1[:, 8], c='darkblue', label='No-Meta')
axs[1].plot(output_array2[:, 0], output_array2[:, 8], c='darkred', label='Meta')
axs[1].set(xlabel='dv_perc [-]', ylabel='Outlet Speed [m/s]')
axs[1].legend(loc='lower right')

axs[2].plot(output_array1[:, 0], output_array1[:, 5], c='darkblue', label='No-Meta')
axs[2].plot(output_array2[:, 0], output_array2[:, 5], c='darkred', label='Meta')
axs[2].set(xlabel='dv_perc [-]', ylabel='Mass Flow Rate [kg/s]')
axs[2].legend(loc='center right')

plt.show()

# %%------------   PLOT RESULTS                         -------------------------------------------------------------> #
fig2, ax2 = plt.subplots()

param = 4

ax2.plot(rotor_array[:, 1], rotor_array[:, param], c='darkblue', label='No-Meta')
ax2.plot(rotor_array2[:, 1], rotor_array2[:, param], c='darkred', label='Meta')

plt.legend()
plt.grid()

plt.show()


# %%------------   PLOT RESULTS                         -------------------------------------------------------------> #
fig2, ax2 = plt.subplots()

param = 8

ax2.plot(rotor_array2[:, 5], rotor_array2[:, param], c='darkred', label='Meta')

plt.legend()
plt.grid()

plt.show()

print('---------Total Enthalpies---------')
print(tesla_turbine2.points[0].get_variable('h'))
print(tesla_turbine2.points[1].get_variable('h'))
print(tesla_turbine2.points[2].get_variable('h'))
print(tesla_turbine2.points[3].get_variable('h'))
print(' ')
print('---------Static Enthalpies---------')
print(tesla_turbine2.static_points[0].get_variable('h'))
print(tesla_turbine2.static_points[1].get_variable('h'))
print(tesla_turbine2.static_points[2].get_variable('h'))
print(tesla_turbine2.static_points[3].get_variable('h'))
print(' ')
print('---------Total Entropies---------')
print(tesla_turbine2.points[0].get_variable('s'))
print(tesla_turbine2.points[1].get_variable('s'))
print(tesla_turbine2.points[2].get_variable('s'))
print(tesla_turbine2.points[3].get_variable('s'))
print(' ')
print('---------Static Entropies---------')
print(tesla_turbine2.static_points[0].get_variable('s'))
print(tesla_turbine2.static_points[1].get_variable('s'))
print(tesla_turbine2.static_points[2].get_variable('s'))
print(tesla_turbine2.static_points[3].get_variable('s'))
print(' ')
print('---------Isentropic Outlet---------')
print(tesla_turbine2.iso_entropic_outlet.get_variable('s'))
print(tesla_turbine2.iso_entropic_outlet.get_variable('h'))

