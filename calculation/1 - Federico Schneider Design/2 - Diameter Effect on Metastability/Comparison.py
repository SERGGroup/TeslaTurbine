# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.multi_phase import TPRotor, TPTeslaGeometry, TPTeslaOptions, TPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

curr_geometry1 = TPTeslaGeometry()
curr_geometry1.stator.d_int = 0.5            # [m]
curr_geometry1.throat_width = 0.003          # [m]
curr_geometry1.rotor.n_channels = 1          # [-]
curr_geometry1.rotor.b_channel = 0.0029      # [m]
curr_geometry1.rotor.roughness = 0.0000005   # [m]

curr_geometry2 = TPTeslaGeometry()
curr_geometry2.stator.d_int = 0.4            # [m]
curr_geometry2.throat_width = 0.003          # [m]
curr_geometry2.rotor.n_channels = 1          # [-]
curr_geometry2.rotor.b_channel = 0.0029      # [m]
curr_geometry2.rotor.roughness = 0.0000005   # [m]

curr_geometry3 = TPTeslaGeometry()
curr_geometry3.stator.d_int = 0.2            # [m]
curr_geometry3.throat_width = 0.003          # [m]
curr_geometry3.rotor.n_channels = 1          # [-]
curr_geometry3.rotor.b_channel = 0.0029      # [m]
curr_geometry3.rotor.roughness = 0.0000005   # [m]

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

tesla_turbine1 = BaseTeslaTurbine("R1234ze", curr_geometry1, curr_options1, stator=TPStatorMil, rotor=TPRotor)
tesla_turbine2 = BaseTeslaTurbine("R1234ze", curr_geometry1, curr_options2, stator=TPStatorMil, rotor=TPRotor)
tesla_turbine3 = BaseTeslaTurbine("R1234ze", curr_geometry2, curr_options1, stator=TPStatorMil, rotor=TPRotor)
tesla_turbine4 = BaseTeslaTurbine("R1234ze", curr_geometry2, curr_options2, stator=TPStatorMil, rotor=TPRotor)
tesla_turbine5 = BaseTeslaTurbine("R1234ze", curr_geometry3, curr_options1, stator=TPStatorMil, rotor=TPRotor)
tesla_turbine6 = BaseTeslaTurbine("R1234ze", curr_geometry3, curr_options2, stator=TPStatorMil, rotor=TPRotor)

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

tesla_turbine3.points[0].set_variable("P", P_in)
tesla_turbine3.points[0].set_variable("x", x_in)
tesla_turbine3.stator.stator_eff = 0.81
tesla_turbine3.rotor.gap_losses_control = True

tesla_turbine4.points[0].set_variable("P", P_in)
tesla_turbine4.points[0].set_variable("x", x_in)
tesla_turbine4.stator.stator_eff = 0.81
tesla_turbine4.rotor.gap_losses_control = True

tesla_turbine5.points[0].set_variable("P", P_in)
tesla_turbine5.points[0].set_variable("x", x_in)
tesla_turbine5.stator.stator_eff = 0.81
tesla_turbine5.rotor.gap_losses_control = True

tesla_turbine6.points[0].set_variable("P", P_in)
tesla_turbine6.points[0].set_variable("x", x_in)
tesla_turbine6.stator.stator_eff = 0.81
tesla_turbine6.rotor.gap_losses_control = True

dv_perc = np.linspace(-0.8, 0.3, 50)

output_array1 = np.empty((len(dv_perc), 9))
output_array2 = np.empty((len(dv_perc), 9))
output_array3 = np.empty((len(dv_perc), 9))
output_array4 = np.empty((len(dv_perc), 9))
output_array5 = np.empty((len(dv_perc), 9))
output_array6 = np.empty((len(dv_perc), 9))

# %%------------   TURBINE-NO_META SOLVING                         --------------------------------------------------> #

for i in tqdm(range(len(dv_perc))):

    tesla_turbine1.rotor.dv_perc = dv_perc[i]

    tesla_turbine1.P_in = P_in
    tesla_turbine1.P_out = P_out
    tesla_turbine1.iterate_pressure()
    rotor_array1 = tesla_turbine1.rotor.get_rotor_array()
    tesla_turbine1.evaluate_performances()

    output_array1[i, 0] = dv_perc[i]
    output_array1[i, 1] = tesla_turbine1.eta_tt
    output_array1[i, 2] = tesla_turbine1.work
    output_array1[i, 3] = tesla_turbine1.power
    output_array1[i, 4] = tesla_turbine1.rotor.rpm
    output_array1[i, 5] = tesla_turbine1.stator.m_dot_s
    output_array1[i, 6] = tesla_turbine1.points[1].get_variable("rho")
    output_array1[i, 7] = tesla_turbine1.static_points[1].get_variable("p")
    output_array1[i, 8] = tesla_turbine1.stator.speed_out.v

    tesla_turbine2.rotor.dv_perc = dv_perc[i]

    tesla_turbine2.P_in = P_in
    tesla_turbine2.P_out = P_out
    tesla_turbine2.iterate_pressure()
    rotor_array2 = tesla_turbine2.rotor.get_rotor_array()
    tesla_turbine2.evaluate_performances()

    output_array2[i, 0] = dv_perc[i]
    output_array2[i, 1] = tesla_turbine2.eta_tt
    output_array2[i, 2] = tesla_turbine2.work
    output_array2[i, 3] = tesla_turbine2.power
    output_array2[i, 4] = tesla_turbine2.rotor.rpm
    output_array2[i, 5] = tesla_turbine2.stator.m_dot_s
    output_array2[i, 6] = tesla_turbine2.points[1].get_variable("rho")
    output_array2[i, 7] = tesla_turbine2.static_points[1].get_variable("p")
    output_array2[i, 8] = tesla_turbine2.stator.speed_out.v


    tesla_turbine3.rotor.dv_perc = dv_perc[i]

    tesla_turbine3.P_in = P_in
    tesla_turbine3.P_out = P_out
    tesla_turbine3.iterate_pressure()
    rotor_array3 = tesla_turbine3.rotor.get_rotor_array()
    tesla_turbine3.evaluate_performances()

    output_array3[i, 0] = dv_perc[i]
    output_array3[i, 1] = tesla_turbine3.eta_tt
    output_array3[i, 2] = tesla_turbine3.work
    output_array3[i, 3] = tesla_turbine3.power
    output_array3[i, 4] = tesla_turbine3.rotor.rpm
    output_array3[i, 5] = tesla_turbine3.stator.m_dot_s
    output_array3[i, 6] = tesla_turbine3.points[1].get_variable("rho")
    output_array3[i, 7] = tesla_turbine3.static_points[1].get_variable("p")
    output_array3[i, 8] = tesla_turbine3.stator.speed_out.v

    tesla_turbine4.rotor.dv_perc = dv_perc[i]

    tesla_turbine4.P_in = P_in
    tesla_turbine4.P_out = P_out
    tesla_turbine4.iterate_pressure()
    rotor_array4 = tesla_turbine4.rotor.get_rotor_array()
    tesla_turbine4.evaluate_performances()

    output_array4[i, 0] = dv_perc[i]
    output_array4[i, 1] = tesla_turbine4.eta_tt
    output_array4[i, 2] = tesla_turbine4.work
    output_array4[i, 3] = tesla_turbine4.power
    output_array4[i, 4] = tesla_turbine4.rotor.rpm
    output_array4[i, 5] = tesla_turbine4.stator.m_dot_s
    output_array4[i, 6] = tesla_turbine4.points[1].get_variable("rho")
    output_array4[i, 7] = tesla_turbine4.static_points[1].get_variable("p")
    output_array4[i, 8] = tesla_turbine4.stator.speed_out.v

    tesla_turbine5.rotor.dv_perc = dv_perc[i]

    tesla_turbine5.P_in = P_in
    tesla_turbine5.P_out = P_out
    tesla_turbine5.iterate_pressure()
    rotor_array5 = tesla_turbine5.rotor.get_rotor_array()
    tesla_turbine5.evaluate_performances()

    output_array5[i, 0] = dv_perc[i]
    output_array5[i, 1] = tesla_turbine5.eta_tt
    output_array5[i, 2] = tesla_turbine5.work
    output_array5[i, 3] = tesla_turbine5.power
    output_array5[i, 4] = tesla_turbine5.rotor.rpm
    output_array5[i, 5] = tesla_turbine5.stator.m_dot_s
    output_array5[i, 6] = tesla_turbine5.points[1].get_variable("rho")
    output_array5[i, 7] = tesla_turbine5.static_points[1].get_variable("p")
    output_array5[i, 8] = tesla_turbine5.stator.speed_out.v

    tesla_turbine6.rotor.dv_perc = dv_perc[i]

    tesla_turbine6.P_in = P_in
    tesla_turbine6.P_out = P_out
    tesla_turbine6.iterate_pressure()
    rotor_array6 = tesla_turbine6.rotor.get_rotor_array()
    tesla_turbine6.evaluate_performances()

    output_array6[i, 0] = dv_perc[i]
    output_array6[i, 1] = tesla_turbine6.eta_tt
    output_array6[i, 2] = tesla_turbine6.work
    output_array6[i, 3] = tesla_turbine6.power
    output_array6[i, 4] = tesla_turbine6.rotor.rpm
    output_array6[i, 5] = tesla_turbine6.stator.m_dot_s
    output_array6[i, 6] = tesla_turbine6.points[1].get_variable("rho")
    output_array6[i, 7] = tesla_turbine6.static_points[1].get_variable("p")
    output_array6[i, 8] = tesla_turbine6.stator.speed_out.v

# %%------------   PLOT RESULTS (1)                     -------------------------------------------------------------> #

# This section plots power and efficiency vs dv_perc for different diameters using metastability

fig, axs = plt.subplots(1, 2, constrained_layout=True)

axs[0].plot(output_array2[:, 0], output_array2[:, 1], c='darkred', label='D = 0.5 m')
axs[0].plot(output_array4[:, 0], output_array4[:, 1], c='darkblue', label='D = 0.4 m')
axs[0].plot(output_array6[:, 0], output_array6[:, 1], c='darkgreen', label='D = 0.2 m')

axs[1].plot(output_array2[:, 0], output_array2[:, 3], c='darkred', label='D = 0.5 m')
axs[1].plot(output_array4[:, 0], output_array4[:, 3], c='darkblue', label='D = 0.4 m')
axs[1].plot(output_array6[:, 0], output_array6[:, 3], c='darkgreen', label='D = 0.2 m')

axs[0].set(xlabel='dv_perc [-]', ylabel='Efficiency [-]', title='Metastability Efficiency')
axs[1].set(xlabel='dv_perc [-]', ylabel='Power [W/ch]', title='Metastability Power')

axs[0].legend()
axs[1].legend()

axs[0].grid()
axs[1].grid()

plt.show()

# %%------------   PLOT RESULTS (2)                     -------------------------------------------------------------> #

# This section plots power and efficiency vs dv_perc for different diameters without using metastability

fig1, axs1 = plt.subplots(1, 2, constrained_layout=True)

axs1[0].plot(output_array1[:, 0], output_array1[:, 1], c='darkred', label='D = 0.5 m')
axs1[0].plot(output_array3[:, 0], output_array3[:, 1], c='darkblue', label='D = 0.4 m')
axs1[0].plot(output_array5[:, 0], output_array5[:, 1], c='darkgreen', label='D = 0.2 m')

axs1[1].plot(output_array1[:, 0], output_array1[:, 3], c='darkred', label='D = 0.5 m')
axs1[1].plot(output_array3[:, 0], output_array3[:, 3], c='darkblue', label='D = 0.4 m')
axs1[1].plot(output_array5[:, 0], output_array5[:, 3], c='darkgreen', label='D = 0.2 m')

axs1[0].set(xlabel='dv_perc [-]', ylabel='Efficiency [-]', title='No-Metastability Efficiency')
axs1[1].set(xlabel='dv_perc [-]', ylabel='Power [W/ch]', title='No-Metastability Power')

axs1[0].legend()
axs1[1].legend()

axs1[0].grid()
axs1[1].grid()

plt.show()

# %%------------   PLOT RESULTS (3)                     -------------------------------------------------------------> #

# This section plots rotor velocity trends for different diameters using metastability

fig2, axs2 = plt.subplots(1, 3, constrained_layout=True)

axs2[0].plot(rotor_array2[:, 1], rotor_array2[:, 2], c='darkred', label='u')
axs2[0].plot(rotor_array2[:, 1], rotor_array2[:, 3], c='darkblue', label='vt')
axs2[0].plot(rotor_array2[:, 1], rotor_array2[:, 4], c='darkgreen', label='vr')

axs2[1].plot(rotor_array4[:, 1], rotor_array4[:, 2], c='darkred', label='u')
axs2[1].plot(rotor_array4[:, 1], rotor_array4[:, 3], c='darkblue', label='vt')
axs2[1].plot(rotor_array4[:, 1], rotor_array4[:, 4], c='darkgreen', label='vr')

axs2[2].plot(rotor_array6[:, 1], rotor_array6[:, 2], c='darkred', label='u')
axs2[2].plot(rotor_array6[:, 1], rotor_array6[:, 3], c='darkblue', label='vt')
axs2[2].plot(rotor_array6[:, 1], rotor_array6[:, 4], c='darkgreen', label='vr')

axs2[0].set(xlabel='Radial DIst [m]', ylabel='Velocity [m/s]', title='D = 0.5 m', ylim=(0, 35))
axs2[1].set(xlabel='Radial DIst [m]', ylabel='Velocity [m/s]', title='D = 0.4 m', ylim=(0, 35))
axs2[2].set(xlabel='Radial DIst [m]', ylabel='Velocity [m/s]', title='D = 0.2 m', ylim=(0, 35))


axs2[0].legend()
axs2[1].legend()
axs2[2].legend()

axs2[0].grid()
axs2[1].grid()
axs2[2].grid()

plt.show()

# %%------------   PLOT RESULTS (4)                     -------------------------------------------------------------> #

# This section plots power and efficiency vs rotational speed for different diameters using metastability

fig3, axs3 = plt.subplots(1, 2, constrained_layout=True)

axs3[0].plot(output_array2[:, 4], output_array2[:, 1], c='darkred', label='D = 0.5 m')
axs3[0].plot(output_array4[:, 4], output_array4[:, 1], c='darkblue', label='D = 0.4 m')
axs3[0].plot(output_array6[:, 4], output_array6[:, 1], c='darkgreen', label='D = 0.2 m')

axs3[1].plot(output_array2[:, 4], output_array2[:, 3], c='darkred', label='D = 0.5 m')
axs3[1].plot(output_array4[:, 4], output_array4[:, 3], c='darkblue', label='D = 0.4 m')
axs3[1].plot(output_array6[:, 4], output_array6[:, 3], c='darkgreen', label='D = 0.2 m')

axs3[0].set(xlabel='RPM [rpm]', ylabel='Efficiency [-]', title='Metastability Efficiency')
axs3[1].set(xlabel='RPM [rp,]', ylabel='Power [W/ch]', title='Metastability Power')

axs3[0].legend()
axs3[1].legend()

axs3[0].grid()
axs3[1].grid()

plt.show()

# %%------------   PLOT RESULTS (5)                     -------------------------------------------------------------> #

# This section plots power and efficiency vs rotational speed for different diameters without using metastability

fig4, axs4 = plt.subplots(1, 2, constrained_layout=True)

axs4[0].plot(output_array1[:, 4], output_array1[:, 1], c='darkred', label='D = 0.5 m')
axs4[0].plot(output_array3[:, 4], output_array3[:, 1], c='darkblue', label='D = 0.4 m')
axs4[0].plot(output_array5[:, 4], output_array5[:, 1], c='darkgreen', label='D = 0.2 m')

axs4[1].plot(output_array1[:, 4], output_array1[:, 3], c='darkred', label='D = 0.5 m')
axs4[1].plot(output_array3[:, 4], output_array3[:, 3], c='darkblue', label='D = 0.4 m')
axs4[1].plot(output_array5[:, 4], output_array5[:, 3], c='darkgreen', label='D = 0.2 m')

axs4[0].set(xlabel='RPM [rpm]', ylabel='Efficiency [-]', title='No-Metastability Efficiency')
axs4[1].set(xlabel='RPM [rpm]', ylabel='Power [W/ch]', title='No-Metastability Power')

axs4[0].legend()
axs4[1].legend()

axs4[0].grid()
axs4[1].grid()

plt.show()

