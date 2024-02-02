# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code import TeslaTurbine, TeslaOptions, TeslaGeometry


# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
curr_geometry = TeslaGeometry()
curr_options = TeslaOptions()
curr_options.stator.iterate_phi = False
curr_geometry.stator.d_int = 0.5    # [m]
P_in = 1000000  # [Pa]
T_in = 400      # K
P_out = 400000  # [Pa]

tt = TeslaTurbine("Air", curr_geometry, curr_options)
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)

tt.static_points[1].set_variable("P", P_out)
tt.static_points[1].set_variable("T", 300)

# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
tt.stator.solve()

# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
tt.rotor.solve()
