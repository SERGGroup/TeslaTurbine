# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.multi_phase import TPRotor, TPTeslaGeometry, TPTeslaOptions, TPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from scipy.interpolate import griddata

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

n_setpoints = 8
dv_perc_list = np.linspace(-0.5, 0.3, n_setpoints)
TW_ref = 0.0015  # [m]
TW_list = np.linspace(0.8, 3, n_setpoints) * TW_ref
dv_perc, TW = np.meshgrid(dv_perc_list, TW_list)

# Input Constraints Data
P_in = 11900000  # [Pa]
P_out = 4500000         # [Pa]
T_in_c = 53.29  # [°C]
T_in = T_in_c + 273.15  # [K]

m_rif = 2            # [kg/s]

output_array_list = list()
rotor_array_max_efficiency_list = list()
rotor_array_max_power_list = list()

Eta = np.empty(len(dv_perc_list) * len(TW_list))
Power = np.empty(len(dv_perc_list) * len(TW_list))

Eta_arr = np.empty(len(dv_perc_list) * len(TW_list))
Power_arr = np.empty(len(dv_perc_list) * len(TW_list))
RPM_arr = np.empty(len(dv_perc_list) * len(TW_list))

output_array = np.empty((len(range(n_setpoints)) * len(range(n_setpoints)), 13))

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
        curr_geometry = TPTeslaGeometry()
        curr_options = TPTeslaOptions()
        curr_options.rotor.integr_variable = 0.04  # [-]
        curr_options.stator.metastability_check = True
        curr_options.rotor.sp_check = False
        curr_options.rotor.n_rotor = 4000

        # Main design Parameters
        curr_geometry.rotor.b_channel = 0.0001
        curr_geometry.rotor.d_ratio = 3  # [m]
        curr_geometry.d_main = 0.10  # [m]

        tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=TPStatorMil, rotor=TPRotor)

        tt.rotor.gap_losses_control = True
        tt.rotor.dv_perc = dv_perc[i,j]

        Z_stat = 3
        tt.geometry.stator.Z_stat = Z_stat

        throat_width = TW[i,j]
        tt.geometry.throat_width = throat_width

        tt.geometry.disc_thickness = 0.0001
        tt.geometry.rotor.n_discs = 1
        tt.geometry.H_s = tt.geometry.rotor.n_discs * (tt.geometry.disc_thickness + tt.geometry.rotor.b_channel)
        tt.m_dot_tot = m_rif
        tt.geometry.n_channels = tt.n_packs

        alpha_in = 85  # [°]
        tt.geometry.alpha1 = alpha_in

        tt.points[0].set_variable("P", P_in)
        tt.points[0].set_variable("T", T_in)
        tt.stator.stator_eff = 0.9
        tt.rotor.gap_losses_control = True

        tt.P_in = P_in
        tt.P_out = P_out
        tt.T_in = T_in

        tt.iterate_pressure()
        tt.evaluate_performances()
        rotor_array = tt.rotor.get_rotor_array()

        if tt.Eta_tesla_ss < 1 and tt.Eta_tesla_ss > 0:

            Eta[j + j_0] = tt.Eta_tesla_ss
            Power[j + j_0] = tt.power * tt.n_packs

            output_array[j + j_0, 0] = tt.rotor.rpm
            output_array[j + j_0, 1] = tt.Eta_tesla_ts
            output_array[j + j_0, 2] = tt.work
            output_array[j + j_0, 3] = tt.power
            output_array[j + j_0, 4] = tt.work2
            output_array[j + j_0, 5] = tt.stator.m_dot_s
            output_array[j + j_0, 6] = tt.n_packs
            output_array[j + j_0, 7] = tt.points[1].get_variable("p")
            output_array[j + j_0, 8] = tt.stator.out_speed
            output_array[j + j_0, 9] = (tt.static_points[2].get_variable("h") - tt.static_points[3].get_variable(
                "h")) / (tt.static_points[0].get_variable("h") - tt.static_points[3].get_variable("h"))
            output_array[j + j_0, 10] = tt.Eta_tesla_ss
            output_array[j + j_0, 11] = tt.points[2].get_variable("p")
            output_array[j + j_0, 12] = tt.volume

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
            output_array[j + j_0, 12] = np.nan

    if tt.Eta_tesla_ss < 1 and tt.Eta_tesla_ss > 0:

        Eta_arr[i * n_setpoints + j] = tt.Eta_tesla_tt2
        Power_arr[i * n_setpoints + j] = tt.power * tt.n_packs
        RPM_arr[i * n_setpoints + j] = tt.rotor.rpm

    else:

        Eta_arr[i * n_setpoints + j] = np.nan
        Power_arr[i * n_setpoints + j] = np.nan
        RPM_arr[i * n_setpoints + j] = np.nan

Eta_res = Eta_arr.reshape(dv_perc.shape)
Power_res = Power_arr.reshape(dv_perc.shape)
RPM_res = RPM_arr.reshape(dv_perc.shape)

 # %%------------   PLOT RESULTS                            ----------------------------------------------------------> #

res1 = (output_array[:, 1]).reshape(dv_perc.shape)
res2 = output_array[:, 5].reshape(dv_perc.shape)
res3 = ((output_array[:, 3] * output_array[:, 6]) / (output_array[:, 12] * 100000)).reshape(dv_perc.shape)
res4 = (output_array[:, 7] / 1000000).reshape(dv_perc.shape)

sigma = 2
res3 = gaussian_filter(res3, sigma)
res1 = gaussian_filter(res1, sigma)


fig, axs = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

ETA = axs[0, 0].contourf(dv_perc, TW, res1, 15)
MFR = axs[0, 1].contourf(dv_perc, TW, res2, 10)
SPEC_POWER = axs[1, 0].contourf(dv_perc, TW, res3, 15)
N_PACKS = axs[1, 1].contourf(dv_perc, TW, res4, 10)

axs[0, 0].set(xlabel='Tip Speed Velocity Ratio [-]', ylabel='Throat Width [m]', title='Efficiency [-]')
axs[0, 1].set(xlabel='Tip Speed Velocity Ratio [-]', ylabel='Throat Width [m]', title='Mass Flow Rate per Channel [kg]')
axs[1, 0].set(xlabel='Tip Speed Velocity Ratio [-]', ylabel='Throat Width [m]', title='Specific Power [kW]')
axs[1, 1].set(xlabel='Tip Speed Velocity Ratio [-]', ylabel='Throat Width [m]', title='Stator Outlet Static Pressure [MPa]')

CB1 = fig.colorbar(ETA, shrink=0.8, ax=axs[0, 0], aspect=30)
CB2 = fig.colorbar(MFR, shrink=0.8, ax=axs[0, 1], aspect=30)
CB3 = fig.colorbar(SPEC_POWER, shrink=0.8, ax=axs[1, 0], aspect=30)
CB4 = fig.colorbar(N_PACKS, shrink=0.8, ax=axs[1, 1], aspect=30)

plt.savefig('grafico.png', dpi=300, bbox_inches='tight')
plt.show()

# %%------------   PLOT RESULTS                            ----------------------------------------------------------> #
method = 'cubic'

# Flatten delle coordinate della griglia originale
x = dv_perc.flatten()
y = TW.flatten()

# Griglia fitta (ad es. 100x100)
xi = np.linspace(x.min(), x.max(), 100)
yi = np.linspace(y.min(), y.max(), 100)
xi, yi = np.meshgrid(xi, yi)

# Interpola ogni risultato (es. res1, res2, etc.)
res1_interp = griddata((x, y), res1.flatten(), (xi, yi), method=method)
res2_interp = griddata((x, y), res2.flatten(), (xi, yi), method=method)
res3_interp = griddata((x, y), res3.flatten(), (xi, yi), method=method)
res4_interp = griddata((x, y), res4.flatten(), (xi, yi), method=method)

# Ridisegna i grafici con i dati interpolati
fig, axs = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

ETA = axs[0, 0].contourf(xi, yi, res1_interp, 15, cmap='viridis')
MFR = axs[0, 1].contourf(xi, yi, res2_interp, 10, cmap='viridis')
SPEC_POWER = axs[1, 0].contourf(xi, yi, res3_interp, 15, cmap='viridis')
N_PACKS = axs[1, 1].contourf(xi, yi, res4_interp, 10, cmap='viridis')

font_axis = 14
font_title = 16

axs[0, 0].set_xlabel('dv_perc [-]', fontsize=font_axis)
axs[0, 0].set_ylabel('Throat Width [m]', fontsize=font_axis)
axs[0, 0].set_title('Efficiency [-]', fontsize=font_title, fontweight='bold')

axs[0, 1].set_xlabel('dv_perc [-]', fontsize=font_axis)
axs[0, 1].set_ylabel('Throat Width [m]', fontsize=font_axis)
axs[0, 1].set_title('Mass Flow Rate per Channel [kg]', fontsize=font_title, fontweight='bold')

axs[1, 0].set_xlabel('dv_perc [-]', fontsize=font_axis)
axs[1, 0].set_ylabel('Throat Width [m]', fontsize=font_axis)
axs[1, 0].set_title('Specific Power [kW/m3]', fontsize=font_title, fontweight='bold')

axs[1, 1].set_xlabel('dv_perc [-]', fontsize=font_axis)
axs[1, 1].set_ylabel('Throat Width [m]', fontsize=font_axis)
axs[1, 1].set_title('Stator Outlet Static Pressure [MPa]', fontsize=font_title, fontweight='bold')

CB1 = fig.colorbar(ETA, shrink=0.8, ax=axs[0, 0], aspect=30)
CB2 = fig.colorbar(MFR, shrink=0.8, ax=axs[0, 1], aspect=30)
CB3 = fig.colorbar(SPEC_POWER, shrink=0.8, ax=axs[1, 0], aspect=30)
CB4 = fig.colorbar(N_PACKS, shrink=0.8, ax=axs[1, 1], aspect=30)

plt.style.use("seaborn-v0_8-whitegrid")

plt.show()