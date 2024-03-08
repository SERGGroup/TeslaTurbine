# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.single_phase import SPStator, SPRotor, SPTeslaGeometry, SPTeslaOptions
from main_code.base_classes import BaseTeslaTurbine

# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
curr_geometry = SPTeslaGeometry()
curr_options = SPTeslaOptions()
curr_options.stator.iterate_phi = True
curr_geometry.stator.d_int = 0.45        # [m]

P_in = 22000000  # [Pa]
T_in = 423.15      # K
P_out = 20000000   # [Pa]

tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStator, rotor=SPRotor)
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)

tt.static_points[1].set_variable("P", P_out)
tt.static_points[1].set_variable("T", 400)


# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
tt.stator.solve()

# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
tt.rotor.solve()