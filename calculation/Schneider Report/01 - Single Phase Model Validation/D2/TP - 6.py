# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.single_phase import SPStator, SPRotor, SPTeslaGeometry, SPTeslaOptions, SPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import matplotlib.pyplot as plt
import numpy as np

# %%------------       INPUT DATA                         -----------------------------------------------------------> #

# Geometric Parameters
curr_geometry = SPTeslaGeometry()
curr_geometry.d_main = 0.217                        # [m]
curr_geometry.throat_width = 0.001                  # [m]
curr_geometry.disc_thickness = 0.0008               # [m]
curr_geometry.alpha1 = 85                           # [°]
curr_geometry.stator.Z_stat = 4                     # [-]

curr_geometry.rotor.gap = 0.001                     # [m]

curr_geometry.rotor.n_discs = 2                     # [-]
curr_geometry.rotor.b_channel = 0.0001              # [m]
curr_geometry.rotor.n_channels = 60                 # [-]
curr_geometry.rotor.d_ratio = 3.92727               # [-]
curr_geometry.rotor.n_packs = 30                    # [-]

single_hs = (curr_geometry.rotor.n_discs - 1) * curr_geometry.disc_thickness + curr_geometry.rotor.n_discs * curr_geometry.rotor.b_channel
curr_geometry.H_s = single_hs * curr_geometry.rotor.n_packs

# Solving Parameters
curr_options = SPTeslaOptions()
curr_options.rotor.integr_variable = 0.03
curr_options.rotor.n_rotor = 400

# Thermodynamic Parameters
P_in = 473535           # [Pa]
T_in = 73.42 + 273.15   # [K]
P_out = 311123          # [Pa]

#Tesla Turbine Initialization
tt = BaseTeslaTurbine("R1233zd(E)", curr_geometry, curr_options, stator=SPStatorMil, rotor=SPRotor)

tt.P_in = P_in
tt.T_in = T_in
tt.P_out = P_out

tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)
tt.points[3].set_variable("P", P_out)

tt.rotor.rpm = 1500
tt.stator.stator_eff = 0.9
tt.rotor.options.sp_check = True

# %%------------       CALCULATIONS                       -----------------------------------------------------------> #

tt.iterate_total_pressure()
rotor_array = tt.rotor.get_rotor_array()
tt.evaluate_performances()

ss = tt.points[1].get_variable("c")
m_dot = tt.stator.m_dot_s
power = tt.power
rpm = tt.rotor.rpm
torque = power / (rpm * 2 * np.pi / 60)

# Power Losses Calculation
A_in_rotor = tt.geometry.rotor.n_discs * 2 * np.pi * tt.geometry.rotor.r_out * tt.geometry.rotor.b_channel
epsilon = 1 - (4 * tt.stator.A_throat) / A_in_rotor
Re2 = tt.points[2].get_variable("rho") * tt.rotor.rotor_inlet_speed.v * tt.geometry.rotor.b_channel / tt.points[2].get_variable("mu")
Cm = 0.003 / Re2 ** 2

T_out = tt.points[2].get_variable("T") - 273.15
Ma_out = tt.stator.out_speed / tt.points[2].get_variable("c")
# C_exp = - 31.08 - 0.023 - 31.08 * T_out + 100.8 * Ma_out + 0.0199 * T_out * Ma_out - 77.04 * Ma_out ** 2
C_exp = 0.05



n_it = 100
P_real = 310
C_up = 0.15
C_down = 0.

for i in range(n_it):

    C_exp = (C_up + C_down) / 2

    # Windage Losses
    P_wind = C_exp * np.pi * tt.geometry.rotor.r_out * tt.geometry.H_s * epsilon * tt.points[2].get_variable(
        "rho") * tt.rotor.first_speed.u ** 3

    # Blockage Losses
    P_block = C_exp * tt.stator.v1_ss * tt.stator.m_dot_s * (
                tt.geometry.rotor.r_out - tt.geometry.rotor.r_int) / tt.geometry.rotor.d_out * tt.rotor.first_speed.u / epsilon

    # Pumping Losses
    P_pump = 4 * Cm * tt.points[2].get_variable("rho") * tt.geometry.rotor.d_out ** 2 * tt.rotor.first_speed.u ** 3

    P_calc = power * tt.geometry.rotor.n_channels - P_wind - P_block - P_pump

    error = abs((P_calc - P_real) / P_real)

    if error < 0.001:
        break
    elif P_calc < P_real:
        C_up = C_exp
    else:
        C_down = C_exp

