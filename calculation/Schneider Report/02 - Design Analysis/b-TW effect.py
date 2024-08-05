# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.multi_phase import TPRotor, TPTeslaGeometry, TPTeslaOptions, TPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

n_setpoints = 10
TW_ref = 0.003  # [m]
b_ref = 0.001  # [m]
TW_list = np.linspace(1, 10, n_setpoints) * TW_ref
b_channel_list = np.linspace(0.3, 1.5, n_setpoints) * b_ref
b_channel, TW = np.meshgrid(b_channel_list, TW_list)

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

Eta = np.empty(len(b_channel_list) * len(TW_list))
Power = np.empty(len(b_channel_list) * len(TW_list))

Eta_arr = np.empty(len(b_channel_list) * len(TW_list))
Power_arr = np.empty(len(b_channel_list) * len(TW_list))
RPM_arr = np.empty(len(b_channel_list) * len(TW_list))

output_array = np.empty((len(range(n_setpoints)) * len(range(n_setpoints)), 10))

# %%------------   CALCULATION                             ----------------------------------------------------------> #

for i in tqdm(range(n_setpoints)):

    print(' ')
    print(' ')
    print('Starting of {}-th b_channel Calculation'.format(i + 1))

    j_0 = i * n_setpoints

    for j in tqdm(range(n_setpoints)):

        # Geometric Parameters
        curr_geometry = TPTeslaGeometry()
        curr_geometry.d_main = 0.2
        curr_geometry.rotor.b_channel = b_channel[i,j]  # [m]
        curr_geometry.disc_thickness = 0.0008  # [m]
        curr_geometry.alpha1 = 85  # [°]
        curr_geometry.stator.Z_stat = 4  # [-]

        curr_geometry.rotor.gap = 0.001             # [m]

        curr_geometry.throat_width = TW[i, j]
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
        new_turbine.m_dot_tot = m_refr

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

Eta_res = Eta_arr.reshape(b_channel.shape)
Power_res = Power_arr.reshape(b_channel.shape)
RPM_res = RPM_arr.reshape(b_channel.shape)

# %%------------   PLOT RESULTS                            ----------------------------------------------------------> #

res1 = output_array[:, 1].reshape(b_channel.shape)
res2 = output_array[:, 5].reshape(b_channel.shape)
res3 = (output_array[:, 3] * output_array[:, 6] / (output_array[:, 0] * 1000)).reshape(b_channel.shape)
res4 = output_array[:, 7].reshape(b_channel.shape)

fig, axs = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

ETA = axs[0, 0].contourf(b_channel, TW, res1, 10)
MFR = axs[0, 1].contourf(b_channel, TW, res2, 8)
SPEC_POWER = axs[1, 0].contourf(b_channel, TW, res3, 10)
PRESSURE = axs[1, 1].contourf(b_channel, TW, res4, 8)

axs[0, 0].set(xlabel='Disk Gap [m]', ylabel='Throat Width [m]', title='Efficiency [-]')
axs[0, 1].set(xlabel='Disk Gap [m]', ylabel='Throat Width [m]', title='Mass Flow Rate per Channel [kg/(s*ch)]')
axs[1, 0].set(xlabel='Disk Gap [m]', ylabel='Throat Width [m]', title='Specific Power [kW / m3]')
axs[1, 1].set(xlabel='Disk Gap [m]', ylabel='Throat Width [m]', title='Stator Outlet Pressure [Pa]')

CB1 = fig.colorbar(ETA, shrink=0.8, ax=axs[0, 0], aspect=30)
CB2 = fig.colorbar(MFR, shrink=0.8, ax=axs[0, 1], aspect=30)
CB3 = fig.colorbar(SPEC_POWER, shrink=0.8, ax=axs[1, 0], aspect=30)
CB4 = fig.colorbar(PRESSURE, shrink=0.8, ax=axs[1, 1], aspect=30)

plt.suptitle("Performance at D = 0.2 m, dv_perc = - 0.4")

plt.show()

# %%------------   PARETO FRONT                            ----------------------------------------------------------> #

fig1, axs1 = plt.subplots(1, 2, constrained_layout=True, figsize=(12, 6))

ETA1 = axs1[0].contourf(b_channel, TW, res1, 10)
SPEC_POWER = axs1[1].contourf(b_channel, TW, res3, 10)

axs1[0].set(xlabel='Disk Gap [m]', ylabel='Throat Width [m]', title='Efficiency [-]')
axs1[1].set(xlabel='Disk Gap [m]', ylabel='Throat Width [m]', title='Specific Power [kW / m3]')

CB11 = fig1.colorbar(ETA1, shrink=0.8, ax=axs1[0], aspect=30)
CB12 = fig1.colorbar(SPEC_POWER, shrink=0.8, ax=axs1[1], aspect=30)

plt.suptitle("Performance at D = 0.2 m, dv_perc = - 0.4")

plt.show()
