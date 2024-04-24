# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
import matplotlib.pyplot as plt
from main_code.sub_classes.multi_phase import TPStator0D, TPRotor, TPTeslaGeometry, TPTeslaOptions, TPStatorMil
from main_code.base_classes import BaseTeslaTurbine

# %%------------   SETUP DATA                             -----------------------------------------------------------> #

curr_geometry = TPTeslaGeometry()
curr_geometry.stator.d_int = 0.5    # [m]
curr_geometry.stator.H_s = 0.0005   # [m]
# curr_geometry.H_s = 0.0005          # [m]  TODO This needs to be cleaned up
curr_geometry.throat_width = 0.003  # [m]
curr_geometry.rotor.n_channels = 1  # [-]
curr_geometry.rotor.b_channel = 0.0029
curr_geometry.rotor.roughness = 0.0000005

curr_options = TPTeslaOptions()
curr_options.stator.iterate_phi = False
curr_options.rotor.profile_rotor = True
curr_options.rotor.sp_check = False
curr_options.rotor.tp_epsilon_model = "sarti"


P_in = 997233  # [Pa]
x_in = 0        # [-]
P_out = 652161  # [Pa]

tt0 = BaseTeslaTurbine("R1234ze", curr_geometry, curr_options, stator=TPStatorMil, rotor=TPRotor)

tt0.points[0].set_variable("P", P_in)
tt0.points[0].set_variable("x", x_in)
tt0.static_points[1].set_variable("P", P_out)
tt0.P_out = P_out
tt0.P_in = P_in
tt0.stator.stator_eff = 0.81

tt0.rotor.omega = 177.9
tt0.rotor.gap_losses_control = True

# %%------------   STATOR SOLVE                           -----------------------------------------------------------> #
tt0.stator.solve()

# %%------------   ROTOR SOLVE                            -----------------------------------------------------------> #
tt0.rotor.solve()
rotor_array = tt0.rotor.get_rotor_array()
tt0.evaluate_performances()

# %%------------   PLOT ROTOR VELOCITIES                        -----------------------------------------------------> #

fig = plt.subplots()

plt.plot(rotor_array[:, 1], rotor_array[:, 3])
plt.plot(rotor_array[:, 1], rotor_array[:, 2])
plt.plot(rotor_array[:, 1], rotor_array[:, 4])
plt.grid()
plt.legend(["vt", "u", "vr"])

plt.show()
