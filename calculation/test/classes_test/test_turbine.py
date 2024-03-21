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

tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStator, rotor=SPRotor)
P_in = 22000000  # [Pa]
tt.P_in = P_in
T_in = 423.15      # K
tt.T_in = T_in
tt.P_out = 19400000   # [Pa]

tt.geometry.d_main = 0.2
tt.geometry.throat_width = 0.0004934
tt.geometry.H_s = 0.00094
tt.geometry.alpha1 = 85

tt.geometry.rotor.d_ratio = 2.5
tt.geometry.rotor.n_channels = 2

tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)
# %%------------       INPUT DATA
tt.iterate_pressure()
rotor_array = tt.rotor.get_rotor_array()
