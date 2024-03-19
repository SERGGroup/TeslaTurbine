# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.multi_phase.stator import TPStator1D, TPStator0D, BaseStator0D, BaseStator1D, BaseStatorStep
from main_code.sub_classes.multi_phase.rotor import TPRotor, TPRotorStep
from main_code.sub_classes.multi_phase import TPTeslaGeometry, TPTeslaOptions
from main_code.base_classes import BaseTeslaTurbine
from REFPROPConnector import ThermodynamicPoint

# %%------------   SETUP DATA                             -----------------------------------------------------------> #
curr_geometry = TPTeslaGeometry()
curr_geometry.stator.d_int = 0.2   # [m]
curr_geometry.stator.H_s = 0.0006   # [m]

curr_options = TPTeslaOptions()
curr_options.stator.iterate_phi = False
curr_options.rotor.profile_rotor = True
curr_options.rotor.tp_epsilon_model = "chisholm"

P_in = 2400000   # [Pa]
x_in = 0         # [-]

tt = BaseTeslaTurbine("r404a", curr_geometry, curr_options, stator=TPStator1D, rotor=TPRotor)
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("x", x_in)

# %%------------   SOLVE                                 -----------------------------------------------------------> #
point = ThermodynamicPoint(["R125", "R143a", "R134a"], [0.44, 0.52, 0.04])
point.set_variable("P", 2.4)
point.set_variable("x", 0)
density = point.get_variable("rho")

tt.stator.m_dot_s = 0.024

# %%------------   SOLVE                                 -----------------------------------------------------------> #
tt.stator.solve_from_flow_rate()

