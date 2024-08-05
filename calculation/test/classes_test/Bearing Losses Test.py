# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.base_classes.support.bearing_losses import BearingLoss
import numpy as np
from main_code.sub_classes.single_phase import SPStator, SPRotor, SPTeslaGeometry, SPTeslaOptions
from main_code.base_classes import BaseTeslaTurbine
import matplotlib.pyplot as plt

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
tt.geometry.alpha_stat = 85
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

M_turb = tt.power / (tt.rotor.rpm * 2 * np.pi / 60)

# %%------------             BEARING LOSSES                ----------------------------------------------------------> #

BL = BearingLoss()

# Geometric Parameters
BL.d_turb = tt.geometry.stator.d_int
BL.d = (tt.geometry.d_main / tt.geometry.rotor.d_ratio)                             # Internal Bearing Diameter
BL.D = BL.d + 0.01                                                                       # External Bearing Diameter
BL.dm = (BL.d + BL.D) / 2
nc = tt.geometry.rotor.n_channels
BL.l_turb = tt.geometry.rotor.b_channel * nc + tt.geometry.disc_thickness * (nc - 1)

# Operational Parameters
BL.nu = 100                             # Lubricating Oil Viscosity in cst (or mm2/s)
BL.n = tt.rotor.rpm

# Thermodynamic Parameters
BL.m_dot = tt.stator.m_dot_s
BL.p = rotor_array[-1, 11]
BL.v = rotor_array[-1, 9]
BL.A_out = np.pi * (tt.geometry.d_main / tt.geometry.rotor.d_ratio) ** 2 / 4

BL.evaluate_bearing_losses()

# %%------------             PLOT RESULTS                 -----------------------------------------------------------> #
fig, axs = plt.subplots(1, 2, constrained_layout=True)

axs[0].plot(rotor_array[:, 1], rotor_array[:, 2], c='darkred', ls='--', label='u')
axs[0].plot(rotor_array[:, 1], rotor_array[:, 3], c='darkblue', ls='--', label='vt')
axs[0].plot(rotor_array[:, 1], rotor_array[:, 4], c='darkgreen', ls='--', label='vr')
axs[0].plot(rotor_array[:, 1], rotor_array[:, 9], c='black', lw = 2.5, label='v')

axs[0].grid()
axs[0].legend()

axs[1].plot(rotor_array[:, 1], rotor_array[:, 2], c='darkred', ls='--', label='u')
axs[1].plot(rotor_array[:, 1], rotor_array[:, 5], c='darkblue', ls='--', label='wt')
axs[1].plot(rotor_array[:, 1], rotor_array[:, 6], c='darkgreen', ls='--', label='wr')
axs[1].plot(rotor_array[:, 1], rotor_array[:, 10], c='black', lw = 2.5, label='w')

axs[1].grid()
axs[1].legend()

plt.show()
