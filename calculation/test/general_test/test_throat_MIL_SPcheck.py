# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.single_phase import SPStator, SPRotor, SPTeslaGeometry, SPTeslaOptions, SPStatorMil
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
from REFPROPConnector import ThermodynamicPoint as TP
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
import matplotlib.pyplot as plt


# %%------------       DATA                        -----------------------------------------------------------> #
n = 100
P_in = 22000000               # [Pa]
T_in = 423.15                 # [K]
P_out = 4000000              # [Pa]
eta_stat = 0.9409             # [-]
H_s = 0.00094                 # [m]
throat_width = 0.0004934      # [m]
z_stat = 4                    # [-]

base_point = TP(["CarbonDioxide"], [1], unit_system="MASS BASE SI")

input_point = base_point.duplicate()
input_point.set_variable("P", P_in)
input_point.set_variable("T", T_in)
H0 = input_point.get_variable("h")
S0 = input_point.get_variable("s")

isentropic_output_point = base_point.duplicate()
isentropic_output_point.set_variable("P", P_out)
isentropic_output_point.set_variable("s", S0)
H1_is = isentropic_output_point.get_variable("h")

v1ss = np.sqrt(2 * (H0 - H1_is))

output_point = base_point.duplicate()
h_out = H0 - eta_stat * (H0 - H1_is)
output_point.set_variable("P", P_out)
output_point.set_variable("h", h_out)

H1 = output_point.get_variable("h")
S1 = output_point.get_variable("S")


# %%------------       MILAZZO STATOR                         -------------------------------------------------------> #
curr_geometry = SPTeslaGeometry()
curr_options = SPTeslaOptions()

tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStatorMil, rotor=SPRotor)
tt.T_in = T_in
tt.P_in = P_in
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)

tt.stator.stator_eff = eta_stat

# Design Parameters
tt.geometry.throat_width = throat_width
tt.geometry.H_s = H_s


# %%------------   CALCULATIONS                        -----------------------------------------------------------> #
Hh = np.linspace(H0, H1, n)
Ph = np.linspace(P_in, P_out, n)
Ss = np.linspace(S0, S0, n)
rhov_throat = 0.

i_max = 0.
v_sound = 0.

v = np.zeros(n)
rhov = np.zeros(n)

m_dot_arr = np.zeros(n)
m_dot_arr_class = np.zeros(n)

int_points = list()

for i in range(n):
    int_points.append(base_point.duplicate())

for i in range(n):

    int_points[i].set_variable("P", Ph[i])
    int_points[i].set_variable("s", Ss[i])

    h_curr = H0 - eta_stat * (H0 - int_points[i].get_variable("h"))

    int_points[i].set_variable("P", Ph[i])
    int_points[i].set_variable("h", h_curr)

    v[i] = np.sqrt(2 * (H0 - h_curr))
    rhov[i] = v[i] * int_points[i].get_variable("rho")

    if rhov[i] > rhov_throat:

        i_max = i
        v_sound = v[i]
        rhov_throat = rhov[i]

    A_throat = throat_width * H_s
    m_dot_s = z_stat * A_throat * rhov_throat
    h_out = H0 - (v_sound ** 2) / 2

    m_dot_arr[i] = m_dot_s

    tt.solve_with_stator_outlet_pressure(Ph[i])
    m_dot_arr_class[i] = tt.stator.m_dot_s

output_point.set_variable("P", P_out)
output_point.set_variable("h", h_out)

v_out = np.sqrt(2 * (H0 - h_out))
Phi = v_out / v1ss


# %%------------       COMPARISON WITH EES                -----------------------------------------------------------> #
P1 = np.linspace(21900000, 21100000, 8)
P2 = np.linspace(21000000, 10000000, 12)
P = np.concatenate((P1, P2))
m = [0.01508, 0.02134, 0.03014, 0.03366, 0.03682, 0.03971, 0.04239, 0.04489, 0.04723, 0.06549, 0.07843, 0.08834, 0.09613, 0.10221, 0.10686, 0.11016, 0.11224, 0.11316, 0.11292, 0.11151]

# %%------------              PLOT 1                       ----------------------------------------------------------> #

fig = plt.subplots()

plt.plot(Ph, m_dot_arr, color='black', linewidth = '1.5')
plt.plot(Ph, m_dot_arr_class, linewidth = '1.5', linestyle = '--')
plt.plot(P[2:], m[2:], linestyle = '--', color='Darkred')
# plt.grid()

plt.xlabel("Back Pressure [Pa]")
plt.ylabel("Mass Flow Rate [kg/s]")

plt.legend(['New Stator Model', 'Old Stator Model'], loc = 'lower left', edgecolor = 'black', facecolor = 'white')
plt.title("Single Phase Validation")
plt.show()
# %%------------              PLOT 2                       ----------------------------------------------------------> #

fig1 = plt.subplots()

plt.plot(v, rhov)
plt.grid()

plt.ylabel("Velocity [m/s]")
plt.xlabel("Density*Velocity [kg/(m2s]")

plt.show()