# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.multi_phase import TPRotor, TPTeslaGeometry, TPTeslaOptions, TPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

curr_geometry = TPTeslaGeometry()
curr_geometry.d_main = 0.2
curr_geometry.rotor.b_channel = 0.0005  # [m]
curr_geometry.disc_thickness = 0.0008  # [m]
curr_geometry.alpha1 = 85  # [°]
curr_geometry.stator.Z_stat = 4  # [-]

curr_geometry.rotor.gap = 0.001  # [m]

curr_geometry.throat_width = 0.015              # [m]
curr_geometry.rotor.n_discs = 1                 # [-]
curr_geometry.rotor.d_ratio = 4                 # [-]
curr_geometry.rotor.n_packs = 1                 # [-]
curr_geometry.rotor.roughness = 0.0000005       # [m]

single_hs = (curr_geometry.rotor.n_discs - 1) * curr_geometry.disc_thickness + curr_geometry.rotor.n_discs * curr_geometry.rotor.b_channel
curr_geometry.H_s = single_hs * curr_geometry.rotor.n_packs

# Solving Parameters
curr_options = TPTeslaOptions()
curr_options.rotor.integr_variable = 0.03
curr_options.rotor.n_rotor = 400

new_turbine = BaseTeslaTurbine('R1234ze', curr_geometry, curr_options, stator=TPStatorMil, rotor=TPRotor)

n_setpoints = 10

x_in = 0                                        # [-]
P_out = 427308                                  # [Pa]
P_in = np.linspace(997233, 1276651, n_setpoints)

rpm_opt = 4224
rpm = np.linspace(0.83, 1.42, n_setpoints) * rpm_opt

output_array = np.empty(((len(P_in) * len(rpm)), 10))

Eta = np.empty(len(rpm) * len(P_in))
Power = np.empty(len(rpm) * len(P_in))
P_eta = np.empty(len(rpm) * len(P_in))
P_power = np.empty(len(rpm) * len(P_in))

max_eff_arr = np.empty((len(rpm)))
max_power_arr = np.empty((len(rpm)))
P_max_eff_arr = np.empty((len(rpm)))
P_max_power_arr = np.empty((len(rpm)))

# %%------------       CALCULATION                         ----------------------------------------------------------> #
for j in tqdm(range(len(rpm))):

    new_turbine.rotor.rpm = rpm[j]

    print(" ")
    print("Starting {}-th RPM calculation".format(j+1))

    max_eff = 0.
    max_power = 0.
    P_max_power = 0.
    P_max_eff = 0.

    for i in range(len(P_in)):

        new_turbine.points[0].set_variable("P", P_in[i])
        new_turbine.points[0].set_variable("x", x_in)
        new_turbine.n_packs = 198
        new_turbine.stator.stator_eff = 0.9
        new_turbine.rotor.gap_losses_control = True

        new_turbine.options.stator.metastability_check = True
        new_turbine.options.rotor.profile_rotor = True
        new_turbine.options.rotor.sp_check = False
        new_turbine.options.rotor.tp_epsilon_model = "sarti"

        new_turbine.P_in = P_in[i]
        new_turbine.P_out = P_out

        new_turbine.iterate_pressure()
        new_turbine.evaluate_performances()

        Eta[j * n_setpoints + i] = new_turbine.Eta_tesla_ss
        Power[j * n_setpoints + i] = new_turbine.power
        P_eta[j * n_setpoints + i] = new_turbine.P_in
        P_power[j * n_setpoints + i] = new_turbine.P_in

        if rpm[j] > rpm_opt:

            if Eta[j * n_setpoints + i] > max_eff or j == 0:
                max_eff = Eta[j * n_setpoints + i]
                P_max_eff = P_eta[j * n_setpoints + i]
            #
            # if Power[j * n_setpoints + i] > max_power or j == 0:
            #     max_power = Power[j * n_setpoints + i]
            #     P_max_power = P_power[j * n_setpoints + i]

        else:
            max_eff = np.nan
            P_max_eff = np.nan
            max_power = np.nan
            P_max_power = np.nan

        output_array[j * n_setpoints + i, 0] = new_turbine.volume
        output_array[j * n_setpoints + i, 1] = new_turbine.Eta_tesla_ss
        output_array[j * n_setpoints + i, 2] = new_turbine.work
        output_array[j * n_setpoints + i, 3] = new_turbine.power
        output_array[j * n_setpoints + i, 4] = new_turbine.rotor.rpm
        output_array[j * n_setpoints + i, 5] = new_turbine.stator.m_dot_s
        output_array[j * n_setpoints + i, 6] = new_turbine.n_packs * new_turbine.power
        output_array[j * n_setpoints + i, 7] = new_turbine.static_points[1].get_variable("p")
        output_array[j * n_setpoints + i, 8] = new_turbine.stator.out_speed
        output_array[j * n_setpoints + i, 9] = new_turbine.P_in

    max_eff_arr[j] = max_eff
    max_power_arr[j] = max_power
    P_max_eff_arr[j] = P_max_eff
    P_max_power_arr[j] = P_max_power

# %%---------------       CALCULATION (OPT)                ----------------------------------------------------------> #

tt = BaseTeslaTurbine('R1234ze', curr_geometry, curr_options, stator=TPStatorMil, rotor=TPRotor)
output_array2 = np.empty((len(P_in), 10))

for i in tqdm(range(len(P_in))):

    tt.rotor.rpm = rpm_opt

    tt.points[0].set_variable("P", P_in[i])
    tt.points[0].set_variable("x", x_in)
    tt.n_packs = 198
    tt.stator.stator_eff = 0.9
    tt.rotor.gap_losses_control = True

    tt.options.stator.metastability_check = True
    tt.options.rotor.profile_rotor = True
    tt.options.rotor.sp_check = False
    tt.options.rotor.tp_epsilon_model = "sarti"

    tt.P_in = P_in[i]
    tt.P_out = P_out

    tt.iterate_pressure()
    tt.evaluate_performances()

    output_array2[i, 0] = tt.volume
    output_array2[i, 1] = tt.Eta_tesla_ss
    output_array2[i, 2] = tt.work
    output_array2[i, 3] = tt.power
    output_array2[i, 4] = tt.rotor.rpm
    output_array2[i, 5] = tt.stator.m_dot_s
    output_array2[i, 6] = tt.n_packs * tt.power
    output_array2[i, 7] = tt.static_points[1].get_variable("p")
    output_array2[i, 8] = tt.stator.out_speed
    output_array2[i, 9] = tt.P_in

# %%------------       PLOT RESULTS                        ----------------------------------------------------------> #
fig, ax = plt.subplots(1, 2, constrained_layout=True, figsize=(9,6))

colors1 = plt.cm.viridis(np.linspace(0.1, 0.9, n_setpoints))
colors2 = plt.cm.viridis(np.linspace(0.1, 0.9, n_setpoints))

for k in range(n_setpoints):

    ax[0].plot(output_array[(k * n_setpoints) : ((k + 1) * n_setpoints - 1), 9],
               output_array[(k * n_setpoints) : ((k + 1) * n_setpoints) - 1, 1], label="rpm = {:.0f}".format(rpm[k]), color=colors1[k])
    ax[1].plot(output_array[(k * n_setpoints): ((k + 1) * n_setpoints)- 1 , 9],
               output_array[(k * n_setpoints): ((k + 1) * n_setpoints) - 1, 6], label="rpm = {:.0f}".format(rpm[k]), color=colors2[k])

ax[0].plot(output_array2[:-1, 9], output_array2[:-1, 1], color='black', linewidth=2.2, linestyle="--")
ax[1].plot(output_array2[:-1, 9], output_array2[:-1, 6], color='black', linewidth=2.2, linestyle="--")

# ax[0].plot(P_max_eff_arr[4:], max_eff_arr[4:], color='red', linewidth=2.2, linestyle="--", label='RPM variable')

ax[0].set(xlabel='Inlet Pressure [Pa]', ylabel='Efficiency [-]')
ax[1].set(xlabel='Inlet Pressure [Pa]', ylabel='Power [W]')

ax[0].legend()
ax[1].legend()

ax[0].grid()
ax[1].grid()

plt.suptitle("OFF-Design - Ambient Temperature")
plt.show()

# %%------------       DP - M_DOT PLOT                     ----------------------------------------------------------> #

m_ref = np.linspace(3.336, 3.336, n_setpoints)

fig, ax = plt.subplots(constrained_layout=True)

colors1 = plt.cm.viridis(np.linspace(0.1, 0.9, n_setpoints))
colors2 = plt.cm.viridis(np.linspace(0.1, 0.9, n_setpoints))

for k in range(n_setpoints):

    ax.plot((output_array[(k * n_setpoints) : ((k + 1) * n_setpoints - 1), 9] / P_out),
               output_array[(k * n_setpoints) : ((k + 1) * n_setpoints) - 1, 5] * 4, label="rpm = {:.0f}".format(rpm[k]), color=colors1[k])

ax.plot(P_in / P_out, m_ref, label='Design Flow Rate', color='black', linestyle="--", linewidth=1.8)

ax.set(xlabel='Beta [-]', ylabel='Mass Flow Rate [kg/s]')

ax.set_xlim([2.34, 2.9])
ax.set_ylim([2, 4.5])

ax.legend()

ax.grid()


plt.suptitle("OFF-Design - Ambient Temperature")
plt.show()