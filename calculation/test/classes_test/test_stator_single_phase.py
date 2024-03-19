# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
import numpy as np
from main_code.sub_classes.single_phase import SPStator, SPRotor, SPTeslaGeometry, SPTeslaOptions
from main_code.base_classes import BaseTeslaTurbine
import matplotlib.pyplot as plt

# %%------------       INPUT DATA                         -----------------------------------------------------------> #
curr_geometry = SPTeslaGeometry()
curr_options = SPTeslaOptions()
curr_options.stator.iterate_phi = True
curr_geometry.stator.d_int = 0.2        # [m]

P_in = 22000000  # [Pa]
T_in = 423.15      # K
P_out = 21000000   # [Pa]

tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStator, rotor=SPRotor)
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)

tt.static_points[1].set_variable("P", P_out)
tt.static_points[1].set_variable("T", 419.15)

tt.geometry.d_main = 0.2
tt.geometry.throat_width = 0.0004934
tt.geometry.H_s = 0.00094
tt.geometry.alpha1 = 85

tt.geometry.rotor.d_ratio = 2.5
tt.geometry.rotor.n_channels = 2

# %%------------     STATOR SOLVE                        -----------------------------------------------------------> #
tt.stator.solve()

# %%------------      ROTOR SOLVE                        -----------------------------------------------------------> #
tt.rotor.solve()
rotor_array = tt.rotor.get_rotor_array()

# %%------------             PLOT                        -----------------------------------------------------------> #

fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

ax.plot(rotor_array[:, 24] *np.pi / 180, rotor_array[:, 1])
ax.set_rmax(0.1)
ax.set_rticks([0.025, 0.05, 0.075, 0.1])
ax.set_rlabel_position(-22.5)
ax.grid(True)

ax.set_title("theta_rel")

plt.show()