# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.base_classes.support.bearing_losses import BearingLoss
import numpy as np
from main_code.sub_classes.single_phase import SPStator, SPRotor, SPTeslaGeometry, SPTeslaOptions
from main_code.base_classes import BaseTeslaTurbine
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

curr_geometry = SPTeslaGeometry()
curr_options = SPTeslaOptions()
curr_options.stator.iterate_phi = True
curr_geometry.stator.d_int = 0.2        # [m]
curr_options.rotor.integr_variable = 0.03
curr_options.stator.metastability_check = False

tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStator, rotor=SPRotor)

P_in = 22000000  # [Pa]
P_out = 19139156
T_in = 423.15      # K
tt.T_in = T_in
tt.P_in = P_in
tt.P_out = P_out
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)

tt.geometry.d_main = 0.2
tt.geometry.H_s = 0.00094
tt.geometry.alpha1 = 85
tt.rotor.gap_losses_control = True

# Design Parameters
tt.geometry.rotor.d_ratio = 2.5
tt.geometry.rotor.b_channel = 0.00007
tt.geometry.throat_width = 0.0004934
tt.geometry.rotor.n_channels = 2
tt.rotor.rpm = 6000

# %%------------             CALCULATIONS                -----------------------------------------------------------> #

tt.iterate_pressure()
rotor_array = tt.rotor.get_rotor_array()
tt.evaluate_performances()

# %%------------             BEARING LOSSES                ----------------------------------------------------------> #

BL = BearingLoss()

BL.m_dot = tt.stator.m_dot
BL.d_turb = tt.geometry.stator.d_int


