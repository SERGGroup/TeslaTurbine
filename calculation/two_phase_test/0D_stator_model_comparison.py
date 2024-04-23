# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
import matplotlib.pyplot as plt
from main_code.sub_classes.multi_phase import TPStator0D, TPRotor, TPTeslaGeometry, TPTeslaOptions, TPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np

# %%------------   SETUP DATA                             -----------------------------------------------------------> #
curr_geometry = TPTeslaGeometry()
curr_geometry.stator.d_int = 0.5    # [m]
curr_geometry.stator.H_s = 0.0005   # [m]
curr_geometry.H_s = 0.0005          # [m]  TODO This needs to be cleaned up
curr_geometry.throat_width = 0.003  # [m]

curr_options = TPTeslaOptions()
curr_options.stator.iterate_phi = True
curr_options.rotor.profile_rotor = True
curr_options.rotor.tp_epsilon_model = "chisholm"

P_in = 997233  # [Pa]
x_in = 0        # [-]
T_in = 44 + 273.15

OUTLET_P = np.linspace(950000, 500000, 100)

tt0 = BaseTeslaTurbine("R1234ze", curr_geometry, curr_options, stator=TPStator0D, rotor=TPRotor)

tt0.points[0].set_variable("P", P_in)
tt0.points[0].set_variable("x", x_in)


tt1 = BaseTeslaTurbine("R1234ze", curr_geometry, curr_options, stator=TPStatorMil, rotor=TPRotor)

tt1.points[0].set_variable("P", P_in)
tt1.points[0].set_variable("T", tt0.points[0].get_variable("T"))

# tt1.points[0].set_variable("T", T_in)


tt1.P_in = tt0.points[0].get_variable("P")  # TODO This needs to be cleaned up too
tt1.T_in = tt0.points[0].get_variable("T")
tt1.stator.stator_eff = 0.81
j = 6
output_array1 = np.empty((len(OUTLET_P), 2 * j - 1))

# %%------------   COMPARE STATOR MODELS                     -------------------------------------------------------> #
for i in range(len(OUTLET_P)):

    tt0.static_points[1].set_variable("P", OUTLET_P[i])
    tt1.static_points[1].set_variable("P", OUTLET_P[i])

    tt0.stator.solve()
    tt1.stator.solve()

    output_array1[i, 0] = OUTLET_P[i]

    output_array1[i, 1] = tt1.stator.m_dot_s
    output_array1[i, 2] = tt1.stator.stator_eff
    output_array1[i, 3] = tt1.stator.out_speed
    output_array1[i, 4] = tt1.stator.stator_mil_out.get_variable("rho")


    output_array1[i, 5] = tt0.stator.m_dot_s
    output_array1[i, 6] = tt0.stator.eta_stat
    output_array1[i, 7] = tt0.stator.v_out
    output_array1[i, 8] = tt0.stator.static_output_point.get_variable("rho")

    output_array1[i, 9] = tt1.stator.phi_n
    output_array1[i, 10] = tt0.stator.phi_n

# %%------------   PLOT COMPARISON                           -------------------------------------------------------> #
fig = plt.subplots()

plt.plot(output_array1[:, 0], output_array1[:, 3], c='darkblue')
plt.plot(output_array1[:, 0], output_array1[:, 7], c='darkred')

plt.grid()
plt.legend(["MILAZZO Model", "OLD Model"])

plt.show()
