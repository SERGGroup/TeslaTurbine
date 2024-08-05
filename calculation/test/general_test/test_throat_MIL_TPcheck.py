# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from REFPROPConnector import ThermodynamicPoint as TP
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# %%------------       DATA                        -----------------------------------------------------------> #
n = 100
P_in = 997233                 # [Pa]
x_in = 0                      # [-]
P_out = 900000                # [Pa]
eta_stat = 0.9               # [-]
H_s = 0.0005                  # [m]
throat_width = 0.003          # [m]
z_stat = 4                    # [-]

base_point = TP(["R1234ze"], [1], unit_system="MASS BASE SI")

input_point = base_point.duplicate()
input_point.set_variable("P", P_in)
input_point.set_variable("x", x_in)
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
T_out = output_point.get_variable("T")
x_out = output_point.get_variable("x")

H1 = output_point.get_variable("h")
S1 = output_point.get_variable("S")

# %%------------   CALCULATIONS                        -----------------------------------------------------------> #

# Hh = np.linspace(H0, H1, n)
Ph = np.linspace(P_in, P_out, n)
Ss = np.linspace(S0, S0, n)
rhov_throat = 0.

i_max = 0.
v_sound = 0.

v = np.zeros(n)
rhov = np.zeros(n)

m_dot_arr = np.zeros(n)

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

output_point.set_variable("P", P_out)
output_point.set_variable("h", h_out)

v_out = np.sqrt(2 * (H0 - h_out))
Phi = v_out / v1ss

# %%------------       COMPARISON WITH EES                -----------------------------------------------------------> #

# %%------------              PLOT                        -----------------------------------------------------------> #

fig, ax1 = plt.subplots(constrained_layout=True)

ax1.plot(Ph, m_dot_arr, label='Mass Flow Rate', color='Darkblue')
ax1.set_xlabel('Back Pressure [Pa]')
ax1.set_ylabel('Mass Flow Rate [kg/s]')
ax1.set_ylim((0, 0.06))
ax1.set_xlim((900000, 1000000))
ax1.xaxis.grid()
ax1.yaxis.grid()

ax2 = ax1.twinx()
ax2.plot(Ph, rhov, label='rho * v', color='Darkred')
ax2.set_ylabel('rho * v [kg/(m2s)]')
ax2.set_ylim((3000, 12000))
ax2.yaxis.set_major_locator(ticker.MultipleLocator(1500))


fig.legend(loc=(0.12, 0.12), edgecolor='black')


plt.show()


