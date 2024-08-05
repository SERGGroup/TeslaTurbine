# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.multi_phase import TPRotor, TPTeslaGeometry, TPTeslaOptions, TPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

n_setpoints = 30
D_ref = 0.35  # [m]
b_ref = 0.001  # [m]
D_turb_list = np.linspace(0.55, 1.43, n_setpoints) * D_ref
b_channel_list = np.linspace(0.3, 1.5, n_setpoints) * b_ref
D_turb, b_channel = np.meshgrid(D_turb_list, b_channel_list)

curr_options = TPTeslaOptions()
curr_options.rotor.profile_rotor = True
curr_options.rotor.sp_check = False
curr_options.rotor.tp_epsilon_model = "sarti"
# curr_options.stator.metastability_check = True

P_in = 997233  # [Pa]
x_in = 0  # [-]
P_out = 427308  # [Pa]

m_refr = 3.336  # [kg/s]

output_array_list = list()
rotor_array_max_efficiency_list = list()
rotor_array_max_power_list = list()

Eta = np.empty(len(D_turb_list) * len(b_channel_list))
Power = np.empty(len(D_turb_list) * len(b_channel_list))

Eta_arr = np.empty(len(D_turb_list) * len(b_channel_list))
Power_arr = np.empty(len(D_turb_list) * len(b_channel_list))
RPM_arr = np.empty(len(D_turb_list) * len(b_channel_list))

output_array = np.empty((len(range(n_setpoints)) * len(range(n_setpoints)), 10))

# %%------------   CALCULATION                             ----------------------------------------------------------> #

for i in tqdm(range(n_setpoints)):

    print(' ')
    print(' ')
    print('Starting of {}-th Diameter Calculation'.format(i))

    j_0 = i * n_setpoints

    for j in tqdm(range(n_setpoints)):

        # Geometric Parameters
        curr_geometry = TPTeslaGeometry()
        curr_geometry.d_main = D_turb[i, j]
        curr_geometry.throat_width = 0.003  # [m]
        curr_geometry.disc_thickness = 0.0008  # [m]
        curr_geometry.alpha1 = 85  # [°]
        curr_geometry.stator.Z_stat = 4  # [-]

        curr_geometry.rotor.gap = 0.001             # [m]

        curr_geometry.rotor.b_channel = b_channel[i, j]
        curr_geometry.rotor.n_discs = 1             # [-]
        curr_geometry.rotor.d_ratio = 3.5           # [-]
        curr_geometry.rotor.n_packs = 1             # [-]
        curr_geometry.rotor.roughness = 0.0000005   # [m]

        single_hs = ( curr_geometry.rotor.n_discs - 1) * curr_geometry.disc_thickness + curr_geometry.rotor.n_discs * curr_geometry.rotor.b_channel
        curr_geometry.H_s = single_hs * curr_geometry.rotor.n_packs

        # Solving Parameters
        curr_options = TPTeslaOptions()
        curr_options.rotor.integr_variable = 0.03
        curr_options.rotor.n_rotor = 400

        new_turbine = BaseTeslaTurbine('R1234ze', curr_geometry, curr_options, stator=TPStatorMil, rotor=TPRotor)
        new_turbine.m_dot_tot = m_refr

        new_turbine.points[0].set_variable("P", P_in)
        new_turbine.points[0].set_variable("x", x_in)
        new_turbine.stator.stator_eff = 0.9
        new_turbine.rotor.gap_losses_control = True

        new_turbine.options.stator.metastability_check = True
        new_turbine.options.rotor.profile_rotor = True
        new_turbine.options.rotor.sp_check = False
        new_turbine.options.rotor.tp_epsilon_model = "sarti"

        new_turbine.rotor.dv_perc = - 0.4

        new_turbine.P_in = P_in
        new_turbine.P_out = P_out
        new_turbine.iterate_pressure()
        new_turbine.evaluate_performances()

        new_turbine.geometry.rotor.n_channels = np.rint(m_refr / new_turbine.stator.m_dot_s)

        Eta[j + j_0] = new_turbine.Eta_tesla_ss
        Power[j + j_0] = new_turbine.power * new_turbine.geometry.rotor.n_channels

        output_array[j + j_0, 0] = new_turbine.volume
        output_array[j + j_0, 1] = new_turbine.Eta_tesla_ss
        output_array[j + j_0, 2] = new_turbine.work
        output_array[j + j_0, 3] = new_turbine.power
        output_array[j + j_0, 4] = new_turbine.rotor.dv_perc
        output_array[j + j_0, 5] = new_turbine.stator.m_dot_s
        output_array[j + j_0, 6] = new_turbine.n_packs
        output_array[j + j_0, 7] = new_turbine.static_points[1].get_variable("p")
        output_array[j + j_0, 8] = new_turbine.stator.out_speed
        output_array[j + j_0, 9] = new_turbine.rotor.out_speed.wt

        Eta_arr[i * n_setpoints + j] = new_turbine.Eta_tesla_ss
        Power_arr[i * n_setpoints + j] = new_turbine.power * new_turbine.geometry.rotor.n_channels
        RPM_arr[i * n_setpoints + j] = new_turbine.rotor.rpm

Eta_res = Eta_arr.reshape(D_turb.shape)
Power_res = Power_arr.reshape(D_turb.shape)
RPM_res = RPM_arr.reshape(D_turb.shape)

# %%------------   PLOT RESULTS                            ----------------------------------------------------------> #

res1 = output_array[:, 1].reshape(D_turb.shape)
res2 = output_array[:, 3].reshape(D_turb.shape)

fig, axs = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)

ETA = axs[0].contourf(D_turb, b_channel, res1)
POWER = axs[1].contourf(D_turb, b_channel, res2)

axs[0].set(xlabel='Diameter [m]', ylabel='b_channel [m]', title='Efficiency [-]')
axs[1].set(xlabel='Diameter [m]', ylabel='b_channel [m]', title='Power [W/ch]')

CB1 = fig.colorbar(ETA, shrink=0.8, ax=axs[0], aspect=30)
CB2 = fig.colorbar(POWER, shrink=0.8, ax=axs[1], aspect=30)

plt.suptitle("Performance at dv_perc = cost")

plt.show()

