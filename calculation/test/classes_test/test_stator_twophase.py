# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.multi_phase import Stator, TPRotor, TPTeslaGeometry, TPTeslaOptions
from main_code.base_classes import BaseTeslaTurbine


# %%------------   SETUP DATA                             -----------------------------------------------------------> #
curr_geometry = TPTeslaGeometry()
curr_geometry.stator.d_int = 0.5    # [m]

curr_options = TPTeslaOptions()
curr_options.stator.iterate_phi = False
curr_options.rotor.profile_rotor = True
curr_options.rotor.tp_epsilon_model = "chisholm"

P_in = 997086  # [Pa]
x_in = 0        # [-]
P_out = 736759  # [Pa]

tt = BaseTeslaTurbine("R1234ze", curr_geometry, curr_options, stator=Stator, rotor=TPRotor)
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("x", x_in)

tt.static_points[1].set_variable("P", P_out)
tt.static_points[1].set_variable("T", 311.7)


# %%------------   SOLVE STATOR                           -----------------------------------------------------------> #
tt.stator.solve()

# %%------------   SOLVE ROTOR                            -----------------------------------------------------------> #
tt.rotor.solve()
