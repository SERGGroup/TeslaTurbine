# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code import TeslaTurbine, TeslaOptions, TeslaGeometry

# %%------------   SETUP DATA                        -----------------------------------------------------------> #
curr_geometry = TeslaGeometry()
curr_options = TeslaOptions()
curr_options.stator.iterate_phi = False
curr_geometry.stator.d_int = 0.5    # [m]
P_in = 997086  # [Pa]
x_in = 0        # [-]
P_out = 736759  # [Pa]

tt = TeslaTurbine("R1234ze", curr_geometry, curr_options)
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("x", x_in)

tt.static_points[1].set_variable("P", P_out)
tt.static_points[1].set_variable("T", 311.7)

T = tt.points[0].get_variable("T")
print(T)
# %%------------   SOLVE                         -----------------------------------------------------------> #
tt.stator.solve()