# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #

from main_code.sub_classes.multi_phase import TPRotor, TPTeslaGeometry, TPTeslaOptions, TPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
import matplotlib.pyplot as plt
from REFPROPConnector import ThermodynamicPoint as TP
from main_code.base_classes.support import Speed, Position
import pandas as pd
import os.path

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

""" VALIDATION PARAMETERS ARE TAKEN FROM:

"Numerical assessment of a two-phase Tesla turbine: Parametric analysis", PH Niknam, L Talluri ...

"""

curr_geometry = TPTeslaGeometry()
curr_options = TPTeslaOptions()
curr_options.rotor.integr_variable = 0.04           # [-]
curr_options.stator.metastability_check = False
curr_options.rotor.sp_check = False
curr_options.rotor.n_rotor = 3000

disk_gap = 0.0004         # [m]
d_ratio = 5.4             # [m]
d_external = 0.216        # [m]

curr_geometry.rotor.b_channel = disk_gap
curr_geometry.rotor.roughness = 0.0000001
curr_geometry.rotor.d_ratio = d_ratio
curr_geometry.d_main = d_external

tt = BaseTeslaTurbine(

    fluid=["R125", "R134A", "R143A"], comp=[0.44, 0.04, 0.52],
    geometry=curr_geometry, options=curr_options,
    stator=TPStatorMil, rotor=TPRotor

)

# CASE 1: Input Thermodynamic Data
P_in = 2281353          # [Pa]
x_in = 0                # [-]

tt.rotor.gap_losses_control = False
tt.rotor.rpm = 2100.8452488130138

inlet_point = TP(["R125", "R134A", "R143A"], [0.44, 0.04, 0.52], unit_system="MASS BASE SI")
inlet_point.set_variable("p", P_in)
inlet_point.set_variable("x", x_in)

inlet_area = np.pi * d_external * disk_gap
m_dot = 0.0016          # [kg/s]
velocity = 26.858       # [m/s]
v_r = m_dot / (inlet_area * inlet_point.get_variable("rho"))

inlet_pos = Position(tt.geometry.d_main / 2, tt.rotor.omega)
inlet_speed = Speed(inlet_pos)
inlet_speed.init_from_codes("v_r", v_r, "v", velocity)

alpha_rad_in = inlet_speed.alpha
alpha_deg_in = alpha_rad_in / np.pi * 180

tt.stator.m_dot_s = m_dot
tt.geometry.rotor.n_channels = 1

# %%--------------------         CALCULATIONS             ----------------------------------------------------------> #

tt.fix_rotor_inlet_condition_quality(P_in, x_in, alpha_rad_in, velocity)
rotor_array = tt.rotor.get_rotor_array()
tt.evaluate_rotor_performances()

parameters = ['it', 'Radius [m]', 'u [m/s]', 'vt [m/s]', 'vr [m/s]', 'wt [m/s]', 'wr [m/s]', 'alpha [rad]', 'beta [rad]',
              'v [m/s]', 'w [m/s]', 'P [Pa]', 'h [J/kg]', 'T [K]', 's [J/kgK]', 'rho [kg/s]', 'mu [Pa*s]', '-', '-', '-',
              'x [-]', '-', '-', '-', 'theta [°]', 'gamma [°]', '-', 'theta90 [°]', 'theta180 [°]', 'theta270 [°]']
rotor_df = pd.DataFrame(rotor_array, columns=parameters)

# RESULTS_FOLDER = os.path.dirname('C:/Users/iggig/PycharmProjects/TeslaTurbine/calculation/00 - Two-Phase Model Validation/RPM Variation')
# excel_file = os.path.join(RESULTS_FOLDER, "results.xlsx")
# rotor_df.to_excel(excel_file, sheet_name='220 rad_s')

power = tt.power
rpm = tt.rotor.rpm
torque = power / (rpm * 2 * np.pi / 60)

v_ratio = rotor_array[-1, 9] / velocity
p_ratio = P_in / rotor_array[-1, 11]
Delta_T = inlet_point.get_variable('T') - rotor_array[-1, 13]

# %%------------             TRAJECTORY PLOT                        -------------------------------------------------> #
fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

ax.plot(rotor_array[:, 25] * np.pi / 180, rotor_array[:, 1])
ax.set_rmax(0.108)
ax.set_rticks([0.025, 0.05, 0.075, 0.1])
ax.set_rlabel_position(-22.5)
ax.grid(True)

ax.set_title("theta_rel")

plt.show()