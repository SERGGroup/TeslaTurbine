# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
import numpy as np
from main_code.sub_classes.single_phase import SPRotor, SPTeslaGeometry, SPTeslaOptions, SPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import matplotlib.pyplot as plt
from REFPROPConnector import ThermodynamicPoint as TP

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

curr_options = SPTeslaOptions()
curr_options.rotor.integr_variable = 0.04           # [-]
curr_options.stator.metastability_check = False
curr_options.rotor.sp_check = True
curr_options.rotor.n_rotor = 3000

# Input Constraints Data
P_in = 10000000         # [bar]
P_out = 3800000         # [bar]
T_in_c = 94             # [°C]
T_in = T_in_c + 273.15  # [K]

# From Ravinath Preliminary Analysis
alpha_in = 85           # [°]
alpha_rad_in = alpha_in * np.pi / 180
throat_width = 0.0005
nozzle_width = 0.02
width = 0.02                                        # [m]
disc_thickness = 0.0001                             # [m]
a_width = throat_width * nozzle_width
Z_stat = 4

inlet_thermopoint = TP(['CarbonDioxide'], [1], unit_system="MASS BASE SI")
inlet_thermopoint.set_variable("P", P_in)
inlet_thermopoint.set_variable("T", T_in)
inlet_ss = inlet_thermopoint.get_variable("c")

mfr_max = inlet_ss * a_width * inlet_thermopoint.get_variable("rho") * Z_stat

test_points = 10
b_channel_array = np.linspace(0.0003, 0.0008, test_points)
dv_perc_array = np.linspace(-0.5, 0.5, test_points)
dv_perc, b_channel = np.meshgrid(dv_perc_array, b_channel_array)

results_array = np.empty((len(dv_perc_array) * len(b_channel_array), 6))
results_array[:] = np.nan

# %%------------   CALCULATIONS                            ----------------------------------------------------------> #

for i in range(len(b_channel_array)):

    curr_geometry = SPTeslaGeometry()
    curr_geometry.rotor.b_channel = b_channel_array[i]
    n_channels = (width + disc_thickness) / (curr_geometry.rotor.b_channel + disc_thickness)
    curr_geometry.rotor.n_channels = n_channels

    curr_geometry.rotor.d_ratio = 3.75                   # [m]
    curr_geometry.d_main = 0.3                           # [m]

    curr_geometry.stator.Z_stat = Z_stat

    tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStatorMil, rotor=SPRotor)
    tt.rotor.gap_losses_control = False

    for j in range(len(dv_perc_array)):

        print(' ')
        print('Calculating {}-th value of b_channel'.format(i))
        print('Calculating {}-th value of dv_perc'.format(j))

        j_0 = i * len(dv_perc_array)

        v_nozzle = 0.
        tt.rotor.dv_perc = dv_perc_array[j]

        mfr_up = mfr_max
        mfr_down = 0.

        print(' ')
        print('Starting Iteration on  Mass Flow Rate -- Limit Values are: MIN {} , MAX {}'.format(0, mfr_max))

        it = 20

        for k in range(it):

            mfr = (mfr_up + mfr_down) / 2
            tt.stator.m_dot_s = mfr
            v_nozzle = mfr / (a_width * inlet_thermopoint.get_variable("rho") * Z_stat)

            tt.fix_rotor_inlet_condition(P_in, T_in, alpha_rad_in, v_nozzle)
            rotor_array = tt.rotor.get_rotor_array()
            tt.evaluate_rotor_performances()

            outlet_pressure = rotor_array[-1, 11]
            error = abs((outlet_pressure - P_out) / P_out)

            if error < 0.0001:
                break
            elif P_out < outlet_pressure:
                mfr_down = mfr
            else:
                mfr_up = mfr

            if abs(mfr - mfr_max) < 0.01:

                print('WARNING!!! The NOZZLE is CHOKED!')
                print(' ')
                break

            print(mfr)

        if tt.power > 0 and not np.isnan(tt.Eta_tesla_ss):

            results_array[j + j_0, 0] = tt.Eta_rotor_ss
            results_array[j + j_0, 1] = tt.power
            results_array[j + j_0, 2] = v_nozzle / inlet_ss
            results_array[j + j_0, 3] = tt.rotor.rpm
            results_array[j + j_0, 4] = tt.power * tt.geometry.rotor.n_channels
            results_array[j + j_0, 5] = mfr


Eta_res = results_array[:, 0].reshape(dv_perc.shape)
Power_res = results_array[:, 4].reshape(dv_perc.shape)
RPM_res = results_array[:, 3].reshape(dv_perc.shape)
m_dot_res = results_array[:, 5].reshape(dv_perc.shape)

# %%------------   PLOT RESULTS                            ----------------------------------------------------------> #

fig, axs = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)

ETA = axs[0].contourf(dv_perc, b_channel, Eta_res)
POWER = axs[1].contourf(dv_perc, b_channel, Power_res)

axs[0].set(xlabel='dv_perc [-]', ylabel='b_channel [m]', title='Efficiency [-]')
axs[1].set(xlabel='dv_perc [-]', ylabel='b_channel [m]', title='Power [W/ch]')

CB1 = fig.colorbar(ETA, shrink=0.8, ax=axs[0], aspect=30)
CB2 = fig.colorbar(POWER, shrink=0.8, ax=axs[1], aspect=30)

plt.show()