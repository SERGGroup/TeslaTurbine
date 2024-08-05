# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
import numpy as np
from main_code.sub_classes.single_phase import SPRotor, SPTeslaGeometry, SPTeslaOptions, SPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import matplotlib.pyplot as plt
from tqdm import tqdm

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

n_setpoints = 20
dv_perc_list = np.linspace(-0.2, 0.4, n_setpoints)
N_packs = 40        # [-]
N_packs_list = np.linspace(0.4, 1.6, n_setpoints) * N_packs
dv_perc, NP = np.meshgrid(dv_perc_list, N_packs_list)

# Input Constraints Data
P_in = 8000000  # [Pa]
# P_out = 3800000         # [Pa]
T_in_c = 94  # [°C]
T_in = T_in_c + 273.15  # [K]

output_array_list = list()
rotor_array_max_efficiency_list = list()
rotor_array_max_power_list = list()

Eta_rotor = np.empty(len(dv_perc_list) * len(N_packs_list))
Power = np.empty(len(dv_perc_list) * len(N_packs_list))

Eta_rotor_arr = np.empty(len(dv_perc_list) * len(N_packs_list))
Power_arr = np.empty(len(dv_perc_list) * len(N_packs_list))
RPM_arr = np.empty(len(dv_perc_list) * len(N_packs_list))

output_array = np.empty((len(range(n_setpoints)) * len(range(n_setpoints)), 8))

Eta_max_arr = np.empty(len(dv_perc_list))
dv_perc_max_arr = np.empty(len(dv_perc_list))
np_max_arr = np.empty(len(dv_perc_list))

# %%------------   CALCULATION                             ----------------------------------------------------------> #

for i in tqdm(range(n_setpoints)):

    print(' ')
    print(' ')
    print('Starting of {}-th dv_perc Calculation'.format(i + 1))

    j_0 = i * n_setpoints

    max_efficiency = 0.
    TW_max = 0.
    dv_perc_max = 0.

    for j in tqdm(range(n_setpoints)):

        # INITIALIZING TURBINE CONDITIONS
        curr_geometry = SPTeslaGeometry()
        curr_options = SPTeslaOptions()
        curr_options.rotor.integr_variable = 0.04  # [-]
        curr_options.stator.metastability_check = False
        curr_options.rotor.sp_check = True
        curr_geometry.rotor.n_channels = 30
        curr_options.rotor.n_rotor = 4000

        # Main design Parameters
        curr_geometry.rotor.b_channel = 0.00005  # [m]
        curr_geometry.rotor.d_ratio = 3  # [m]
        curr_geometry.d_main = 0.15  # [m]

        tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStatorMil, rotor=SPRotor)

        tt.rotor.gap_losses_control = False
        tt.rotor.dv_perc = dv_perc[i,j]

        Z_stat = 1
        tt.geometry.stator.Z_stat = Z_stat

        throat_width = 0.0035

        tt.geometry.throat_width = throat_width
        tt.geometry.disc_thickness = 0.0001
        tt.geometry.rotor.n_discs = 1
        tt.geometry.H_s = 0.00005
        tt.n_packs = NP[i,j]
        tt.geometry.n_channels = tt.n_packs

        nozzle_width = tt.geometry.H_s * tt.n_packs
        a_width = throat_width * nozzle_width

        inlet_thermopoint = tt.points[0].duplicate()
        inlet_thermopoint.set_variable("P", P_in)
        inlet_thermopoint.set_variable("T", T_in)
        inlet_ss = inlet_thermopoint.get_variable("c")

        mfr_max = inlet_ss * a_width * inlet_thermopoint.get_variable("rho") * Z_stat

        mfr = 0.1  # [kg/s]
        tt.stator.m_dot_s = mfr

        inlet_tp_it = tt.points[0].duplicate()
        inlet_tp_it.set_variable("p", P_in)
        inlet_tp_it.set_variable("h", inlet_thermopoint.get_variable("h"))

        rho_it = inlet_tp_it.get_variable("rho")
        T_it = inlet_tp_it.get_variable("T")

        A = np.pi * tt.geometry.d_main * tt.geometry.rotor.b_channel
        ur = mfr / (rho_it * A)

        alpha_in = 85  # [°]
        alpha_rad_in = alpha_in * np.pi / 180

        v_nozzle = mfr / (rho_it * a_width * Z_stat)

        #SOLVING
        tt.fix_rotor_inlet_condition(P_in, T_in, alpha_rad_in, v_nozzle)
        rotor_array = tt.rotor.get_rotor_array()
        tt.evaluate_rotor_performances()

        #DATA COLLECTING
        Eta_rotor[j + j_0] = tt.Eta_rotor_tt
        Power[j + j_0] = tt.power * tt.geometry.rotor.n_channels

        if Eta_rotor[j + j_0] > max_efficiency or j == 0:
            rotor_array_max_efficiency = tt.rotor.get_rotor_array()
            max_efficiency = Eta_rotor[j + j_0]
            NP_max = NP[i,j]
            dv_perc_max = dv_perc[i,j]

        output_array[j + j_0, 0] = tt.volume
        output_array[j + j_0, 1] = tt.Eta_rotor_tt
        output_array[j + j_0, 2] = rotor_array[-1, 11]
        output_array[j + j_0, 3] = tt.power * tt.n_packs
        output_array[j + j_0, 4] = tt.rotor.rpm
        output_array[j + j_0, 5] = tt.rotor.rotor_inlet_speed.v
        output_array[j + j_0, 6] = (P_in - rotor_array[-1, 11])/100000
        output_array[j + j_0, 7] = v_nozzle

    Eta_rotor_arr[j + j_0] = tt.Eta_rotor_tt
    Power_arr[j + j_0] = tt.power * tt.geometry.rotor.n_channels
    RPM_arr[j + j_0] = tt.rotor.rpm

    Eta_max_arr[i] = max_efficiency
    np_max_arr[i] = NP_max
    dv_perc_max_arr[i] = dv_perc_max

Eta_rotor_res = Eta_rotor.reshape(dv_perc.shape)
Power_res = Power_arr.reshape(dv_perc.shape)
RPM_res = RPM_arr.reshape(dv_perc.shape)

# %%------------   PLOT RESULTS                            ----------------------------------------------------------> #

res1 = output_array[:, 1].reshape(dv_perc.shape)
res2 = output_array[:, 7].reshape(dv_perc.shape)
res3 = output_array[:, 3].reshape(dv_perc.shape)
res4 = output_array[:, 6].reshape(dv_perc.shape)

fig, axs = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

ETA = axs[0, 0].contourf(dv_perc, NP, res1, 15)
RPM = axs[0, 1].contourf(dv_perc, NP, res2, 8)
POWER = axs[1, 0].contourf(dv_perc, NP, res3, 10)
PR = axs[1, 1].contourf(dv_perc, NP, res4, 8)

# axs[0, 0].plot(dv_perc_max_arr, np_max_arr, color='Darkred', linewidth='3', linestyle='-', alpha=0.5)
axs[0, 0].set(xlabel='dv_perc [-]', ylabel='Packs Number [-]', title='Rotor Efficiency [-]')
axs[0, 1].set(xlabel='dv_perc [-]', ylabel='Packs Number [-]', title='Nozzle Outlet Velocity [m/s]')
axs[1, 0].set(xlabel='dv_perc [-]', ylabel='Packs Number [-]', title='Power [W]')
axs[1, 1].set(xlabel='dv_perc [-]', ylabel='Packs Number [-]', title='Static-to-Static Pressure Variation [bar]')

CB1 = fig.colorbar(ETA, shrink=0.8, ax=axs[0, 0], aspect=30)
CB2 = fig.colorbar(RPM, shrink=0.8, ax=axs[0, 1], aspect=30)
CB3 = fig.colorbar(POWER, shrink=0.8, ax=axs[1, 0], aspect=30)
CB4 = fig.colorbar(PR, shrink=0.8, ax=axs[1, 1], aspect=30)

plt.suptitle("Performance at P_in = 80 bar, TW = 3.5 mm, b = 0.05 mm")

plt.show()