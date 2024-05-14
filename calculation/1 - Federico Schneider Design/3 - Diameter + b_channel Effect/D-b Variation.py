# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.multi_phase import TPRotor, TPTeslaGeometry, TPTeslaOptions, TPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt


# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

n_setpoints = 3
D_ref = 0.35     # [m]
b_ref = 0.001   # [m]
D_turb_list = np.linspace(0.55, 1.45, n_setpoints) * D_ref
b_channel_list = np.linspace(0.5, 3, n_setpoints) * b_ref
D_turb, b_channel = np.meshgrid(D_turb_list, b_channel_list)

curr_geometry = TPTeslaGeometry()
curr_geometry.stator.d_int = D_ref          # [m]
curr_geometry.rotor.b_channel = b_ref       # [m]
curr_geometry.throat_width = 0.003          # [m]
curr_geometry.rotor.n_channels = 1          # [-]
curr_geometry.rotor.roughness = 0.0000005   # [m]

curr_options = TPTeslaOptions()
curr_options.rotor.profile_rotor = True
curr_options.rotor.sp_check = False
curr_options.rotor.tp_epsilon_model = "sarti"
curr_options.stator.metastability_check = True

# ref_turbine = BaseTeslaTurbine('R1234ze', curr_geometry, curr_options, stator=TPStatorMil, rotor=TPRotor)
turbine_list = list()

for i in range(n_setpoints):

    for j in range(n_setpoints):

        curr_geometry.stator.d_int = D_turb[i, j]
        curr_geometry.rotor.b_channel = b_channel[i, j]
        new_turbine = BaseTeslaTurbine('R1234ze', curr_geometry, curr_options, stator=TPStatorMil, rotor=TPRotor)

        print("---------New Turbine---------")
        print(new_turbine.geometry.stator.d_int)
        print(new_turbine.geometry.rotor.b_channel)
        print(" ")

        turbine_list.append(new_turbine)

        print("---------List Element---------")
        print(turbine_list[i * n_setpoints + j].geometry.stator.d_int)
        print(turbine_list[i * n_setpoints + j].geometry.rotor.b_channel)
        print(" ")

print(turbine_list[0].geometry.stator.d_int)
print(turbine_list[1].geometry.stator.d_int)
print(turbine_list[2].geometry.stator.d_int)
print(turbine_list[3].geometry.stator.d_int)

P_in = 997233         # [Pa]
x_in = 0              # [-]
P_out = 427308        # [Pa]

dv_perc = np.linspace(-0.8, 0.3, 4)

output_array_list = list()
rotor_array_max_efficiency_list = list()
rotor_array_max_power_list = list()

Eta = np.empty(len(dv_perc) * len(turbine_list))
Power = np.empty(len(dv_perc) * len(turbine_list))

Eta_max_arr = np.empty(len(turbine_list))
Power_max_arr = np.empty(len(turbine_list))


# %%------------   CALCULATION                             ----------------------------------------------------------> #

for k in tqdm(range(len(turbine_list))):

    p_0 = k * len(dv_perc)
    max_efficiency = 0.
    max_power = 0.
    output_array = np.empty((len(dv_perc), 9))

    for p in range(len(dv_perc)):

        turbine_list[k].rotor.dv_perc = dv_perc[p]

        turbine_list[k].points[0].set_variable("P", P_in)
        turbine_list[k].points[0].set_variable("x", x_in)
        turbine_list[k].stator.stator_eff = 0.81
        turbine_list[k].rotor.gap_losses_control = True

        turbine_list[k].P_in = P_in
        turbine_list[k].P_out = P_out
        turbine_list[k].iterate_pressure()
        turbine_list[k].evaluate_performances()

        Eta[p + p_0] = turbine_list[k].eta_tt
        Power[p + p_0] = turbine_list[k].power

        if Eta[p + p_0] > max_efficiency or p == 0:
            rotor_array_max_efficiency = turbine_list[k].rotor.get_rotor_array()
            max_efficiency = Eta[p + p_0]

        if Power[p + p_0] > max_power or p == 0:
            rotor_array_max_power = turbine_list[k].rotor.get_rotor_array()
            max_power = Power[p + p_0]

        output_array[p, 0] = dv_perc[p]
        output_array[p, 1] = turbine_list[k].eta_tt
        output_array[p, 2] = turbine_list[k].work
        output_array[p, 3] = turbine_list[k].power
        output_array[p, 4] = turbine_list[k].rotor.rpm
        output_array[p, 5] = turbine_list[k].stator.m_dot_s
        output_array[p, 6] = turbine_list[k].points[1].get_variable("rho")
        output_array[p, 7] = turbine_list[k].static_points[1].get_variable("p")
        output_array[p, 8] = turbine_list[k].stator.speed_out.v

    output_array_list.append(output_array)
    rotor_array_max_efficiency_list.append(rotor_array_max_efficiency)
    rotor_array_max_power_list.append(rotor_array_max_power)
    Eta_max_arr[k] = max_efficiency
    Power_max_arr[k] = max_power

Eta_res = Eta_max_arr.reshape(D_turb.shape)
Power_res = Power_max_arr.reshape(D_turb.shape)

# %%------------   PLOT RESULTS                            ----------------------------------------------------------> #

fig, axs = plt.subplots(1, 2, constrained_layout=True)

axs[0].contourf(D_turb, b_channel, Eta_res)
axs[1].contourf(D_turb, b_channel, Power_res)

axs[0].set(xlabel='Diameter [m]', ylabel='b_channel [m]', title='Maximum Efficiency [-]')
axs[1].set(xlabel='Diameter [m]', ylabel='b_channel [m]', title='Maximum Power [W/ch]')

plt.show()
