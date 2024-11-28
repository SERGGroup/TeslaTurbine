# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.single_phase import SPRotor, SPTeslaGeometry, SPTeslaOptions, SPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

# The design parameters analyzed here are throat width and dv_perc. This section of code define their range of variation
n_setpoints = 8
dv_perc_list = np.linspace(-0.4, 0.2, n_setpoints)
TW_ref = 0.00114  # [m]
TW_list = np.linspace(0.8, 1.2, n_setpoints) * TW_ref
dv_perc, TW = np.meshgrid(dv_perc_list, TW_list)

# Input Thermodynamic Data
P_in = 10100000  # [Pa]
P_out = 3900000         # [Pa]
T_in_c = 94  # [°C]
T_in = T_in_c + 273.15  # [K]

m_rif = 0.1             # [kg/s]

# Initializing output arrays
output_array_list = list()
rotor_array_max_efficiency_list = list()
rotor_array_max_power_list = list()

Eta = np.empty(len(dv_perc_list) * len(TW_list))
Power = np.empty(len(dv_perc_list) * len(TW_list))

Eta_arr = np.empty(len(dv_perc_list) * len(TW_list))
Power_arr = np.empty(len(dv_perc_list) * len(TW_list))
RPM_arr = np.empty(len(dv_perc_list) * len(TW_list))

output_array = np.empty((len(range(n_setpoints)) * len(range(n_setpoints)), 12))

Eta_max_arr = np.empty(len(dv_perc_list))
D_max_arr = np.empty(len(dv_perc_list))
b_max_arr = np.empty(len(dv_perc_list))


# %%------------   CALCULATION                             ----------------------------------------------------------> #

for i in tqdm(range(n_setpoints)):

    print(' ')
    print(' ')
    print('Starting of {}-th dv_perc Calculation'.format(i + 1))

    j_0 = i * n_setpoints

    max_efficiency = 0.
    D_max = 0.
    b_max = 0.

    for j in tqdm(range(n_setpoints)):

        # INITIALIZING TURBINE CONDITIONS
        curr_geometry = SPTeslaGeometry()
        curr_options = SPTeslaOptions()
        curr_options.rotor.integr_variable = 0.04  # [-]
        curr_options.stator.metastability_check = False
        curr_options.rotor.sp_check = True
        curr_options.rotor.n_rotor = 4000

        # These are the other design parameters. In this case, they're kept constant (only tw and dv_perc are variating)
        curr_geometry.rotor.b_channel = 0.00005 # [m]
        curr_geometry.rotor.d_ratio = 3  # [m]
        curr_geometry.d_main = 0.15  # [m]

        # Initializing the turbine
        tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStatorMil, rotor=SPRotor)

        tt.rotor.dv_perc = dv_perc[i,j]

        Z_stat = 1
        tt.geometry.stator.Z_stat = Z_stat

        throat_width = TW[i,j]
        tt.geometry.throat_width = throat_width

        tt.geometry.disc_thickness = 0.0001
        tt.geometry.rotor.n_discs = 1
        tt.geometry.H_s = 0.00015
        tt.m_dot_tot = m_rif
        tt.geometry.n_channels = tt.n_packs

        alpha_in = 89  # [°]
        tt.geometry.alpha1 = alpha_in

        tt.points[0].set_variable("P", P_in)
        tt.points[0].set_variable("T", T_in)
        tt.stator.stator_eff = 0.9
        tt.rotor.gap_losses_control = True

        tt.P_in = P_in
        tt.P_out = P_out
        tt.T_in = T_in

        # This is the section where the calculation is performed
        tt.iterate_pressure()
        tt.evaluate_performances()
        rotor_array = tt.rotor.get_rotor_array()

        # This is just a check conditions to avoid numerical problems (i.e. negative efficiency or higher than 100%
        if tt.Eta_tesla_ss < 1 and tt.Eta_tesla_ss > 0:

            Eta[j + j_0] = tt.Eta_tesla_ss
            Power[j + j_0] = tt.power * tt.n_packs

            output_array[j + j_0, 0] = tt.rotor.rpm
            output_array[j + j_0, 1] = tt.Eta_tesla_tt2
            output_array[j + j_0, 2] = tt.work
            output_array[j + j_0, 3] = tt.power
            output_array[j + j_0, 4] = tt.work2
            output_array[j + j_0, 5] = tt.stator.m_dot_s
            output_array[j + j_0, 6] = tt.n_packs
            output_array[j + j_0, 7] = tt.static_points[1].get_variable("p")
            output_array[j + j_0, 8] = tt.stator.out_speed
            output_array[j + j_0, 9] = (tt.static_points[2].get_variable("h") - tt.static_points[3].get_variable("h")) / (tt.static_points[0].get_variable("h") - tt.static_points[3].get_variable("h"))
            output_array[j + j_0, 10] = tt.Eta_tesla_ss
            output_array[j + j_0, 11] = tt.Eta_tesla_ts

        else:
            output_array[j + j_0, 0] = np.nan
            output_array[j + j_0, 1] = np.nan
            output_array[j + j_0, 2] = np.nan
            output_array[j + j_0, 3] = np.nan
            output_array[j + j_0, 4] = np.nan
            output_array[j + j_0, 5] = np.nan
            output_array[j + j_0, 6] = np.nan
            output_array[j + j_0, 7] = np.nan
            output_array[j + j_0, 8] = np.nan
            output_array[j + j_0, 9] = np.nan
            output_array[j + j_0, 10] = np.nan
            output_array[j + j_0, 11] = np.nan

    # Exporting results in the arrays
    if tt.Eta_tesla_ss < 1 and tt.Eta_tesla_ss > 0:

        Eta_arr[i * n_setpoints + j] = tt.Eta_tesla_tt2
        Power_arr[i * n_setpoints + j] = tt.power * tt.n_packs
        RPM_arr[i * n_setpoints + j] = tt.rotor.rpm

    else:

        Eta_arr[i * n_setpoints + j] = np.nan
        Power_arr[i * n_setpoints + j] = np.nan
        RPM_arr[i * n_setpoints + j] = np.nan

# Reshaping arrays for the final plots
Eta_res = Eta_arr.reshape(dv_perc.shape)
Power_res = Power_arr.reshape(dv_perc.shape)
RPM_res = RPM_arr.reshape(dv_perc.shape)

# %%------------   PLOT RESULTS                            ----------------------------------------------------------> #

res1 = output_array[:, 1].reshape(dv_perc.shape)
res2 = output_array[:, 0].reshape(dv_perc.shape)
res3 = ((output_array[:, 3] * output_array[:, 6]) / 1000).reshape(dv_perc.shape)
res4 = (output_array[:, 6]).reshape(dv_perc.shape)

sigma = 1.5
res3 = gaussian_filter(res3, sigma)
res4 = gaussian_filter(res4, sigma)

fig, axs = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

ETA = axs[0, 0].contourf(dv_perc, TW, res1, 10)
NPACK = axs[0, 1].contourf(dv_perc, TW, res2, 8)
SPEC_POWER = axs[1, 0].contourf(dv_perc, TW, res3, 10)
PRESSURE = axs[1, 1].contourf(dv_perc, TW, res4, 10)

axs[0, 0].set(xlabel='dv_perc [-]', ylabel='Throat Width [m]', title='Efficiency [-]')
axs[0, 1].set(xlabel='dv_perc [-]', ylabel='Throat Width [m]', title='Rotational Speed [rpm]')
axs[1, 0].set(xlabel='dv_perc [-]', ylabel='Throat Width [m]', title='Power [kW]')
axs[1, 1].set(xlabel='dv_perc [-]', ylabel='Throat Width [m]', title='Number of Channels [-]')

CB1 = fig.colorbar(ETA, shrink=0.8, ax=axs[0, 0], aspect=30)
CB2 = fig.colorbar(NPACK, shrink=0.8, ax=axs[0, 1], aspect=30)
CB3 = fig.colorbar(SPEC_POWER, shrink=0.8, ax=axs[1, 0], aspect=30)
CB4 = fig.colorbar(PRESSURE, shrink=0.8, ax=axs[1, 1], aspect=30)

plt.suptitle("Performance at D = 0.15 m, b = 0.05 mm, D_ratio = 3", fontsize='16')

plt.show()

# %%------------             TRAJECTORY PLOT                        -------------------------------------------------> #
fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

ax.plot(rotor_array[:, 25] * np.pi / 180, rotor_array[:, 1])
ax.set_rmax(0.075)
ax.set_rticks([tt.geometry.d_main / tt.geometry.rotor.d_ratio /2,
               (tt.geometry.d_main - tt.geometry.d_main / tt.geometry.rotor.d_ratio) / 2])
ax.grid(True)

ax.set_title("theta_rel")

plt.show()
