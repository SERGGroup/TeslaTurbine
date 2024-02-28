# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code import TeslaTurbine, TeslaOptions, TeslaGeometry


# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
curr_geometry = TeslaGeometry()
curr_options = TeslaOptions()
curr_options.stator.iterate_phi = False
curr_geometry.stator.d_int = 0.45    # [m]
P_in = 22000000  # [Pa]
T_in = 423.15      # K
P_out = 20000000   # [Pa]

tt = TeslaTurbine("Air", curr_geometry, curr_options)
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)

tt.static_points[1].set_variable("P", P_out)
tt.static_points[1].set_variable("T", 300)


# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
tt.stator.solve()
