# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.multi_phase import TPRotor, TPTeslaGeometry, TPTeslaOptions, TPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt


# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

n_setpoints = 5
D_ref = 0.35     # [m]
b_ref = 0.001   # [m]
D_turb_list = np.linspace(0.55, 1.45, n_setpoints) * D_ref
b_channel_list = np.linspace(0.5, 3, n_setpoints) * b_ref
D_turb, b_channel = np.meshgrid(D_turb_list, b_channel_list)

curr_options = TPTeslaOptions()
curr_options.rotor.profile_rotor = True
curr_options.rotor.sp_check = False
curr_options.rotor.tp_epsilon_model = "sarti"
curr_options.stator.metastability_check = False

P_in = 997233         # [Pa]
x_in = 0              # [-]
P_out = 427308        # [Pa]

dv_perc = np.linspace(-0.8, 0.3, 15)

output_array_list = list()
rotor_array_max_efficiency_list = list()
rotor_array_max_power_list = list()

Eta = np.empty(len(dv_perc) * len(D_turb_list) * len(b_channel_list))
Power = np.empty(len(dv_perc) * len(D_turb_list) * len(b_channel_list))

Eta_max_arr = np.empty(len(D_turb_list) * len(b_channel_list))
Power_max_arr = np.empty(len(D_turb_list) * len(b_channel_list))

# %%------------   CALCULATION                             ----------------------------------------------------------> #

for i in tqdm(range(n_setpoints)):

    for j in range(n_setpoints):

        curr_geometry = TPTeslaGeometry()
        curr_geometry.throat_width = 0.003  # [m]
        curr_geometry.rotor.n_channels = 1  # [-]
        curr_geometry.rotor.roughness = 0.0000005  # [m]

        curr_geometry.stator.d_int = D_turb[i, j]
        curr_geometry.rotor.b_channel = b_channel[i, j]

        new_turbine = BaseTeslaTurbine('R1234ze', curr_geometry, curr_options, stator=TPStatorMil, rotor=TPRotor)

        p_0 = (i * n_setpoints + j) * len(dv_perc)
        max_efficiency = 0.
        max_power = 0.
        output_array = np.empty((len(dv_perc), 9))

        for p in range(len(dv_perc)):

            new_turbine.rotor.dv_perc = dv_perc[p]

            new_turbine.points[0].set_variable("P", P_in)
            new_turbine.points[0].set_variable("x", x_in)
            new_turbine.stator.stator_eff = 0.81
            new_turbine.rotor.gap_losses_control = True

            new_turbine.P_in = P_in
            new_turbine.P_out = P_out
            new_turbine.iterate_pressure()
            new_turbine.evaluate_performances()

            Eta[p + p_0] = new_turbine.Eta_tesla_ss
            Power[p + p_0] = new_turbine.power

            if Eta[p + p_0] > max_efficiency or p == 0:
                rotor_array_max_efficiency = new_turbine.rotor.get_rotor_array()
                max_efficiency = Eta[p + p_0]

            if Power[p + p_0] > max_power or p == 0:
                rotor_array_max_power = new_turbine.rotor.get_rotor_array()
                max_power = Power[p + p_0]

            output_array[p, 0] = dv_perc[p]
            output_array[p, 1] = new_turbine.Eta_tesla_ss
            output_array[p, 2] = new_turbine.work
            output_array[p, 3] = new_turbine.power
            output_array[p, 4] = new_turbine.rotor.rpm
            output_array[p, 5] = new_turbine.stator.m_dot_s
            output_array[p, 6] = new_turbine.points[1].get_variable("rho")
            output_array[p, 7] = new_turbine.static_points[1].get_variable("p")
            output_array[p, 8] = new_turbine.stator.out_speed

        output_array_list.append(output_array)
        rotor_array_max_efficiency_list.append(rotor_array_max_efficiency)
        rotor_array_max_power_list.append(rotor_array_max_power)
        Eta_max_arr[i * n_setpoints + j] = max_efficiency
        Power_max_arr[i * n_setpoints + j] = max_power

Eta_res = Eta_max_arr.reshape(D_turb.shape)
Power_res = Power_max_arr.reshape(D_turb.shape)


# %%------------   PLOT RESULTS                            ----------------------------------------------------------> #

fig, axs = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)

ETA = axs[0].contourf(D_turb, b_channel, Eta_res)
POWER = axs[1].contourf(D_turb, b_channel, Power_res)

axs[0].set(xlabel='Diameter [m]', ylabel='b_channel [m]', title='Maximum Efficiency [-]')
axs[1].set(xlabel='Diameter [m]', ylabel='b_channel [m]', title='Maximum Power [W/ch]')

CB1 = fig.colorbar(ETA, shrink=0.8, ax=axs[0], aspect=30)
CB2 = fig.colorbar(POWER, shrink=0.8, ax=axs[1], aspect=30)

plt.show()



