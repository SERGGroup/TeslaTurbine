# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from REFPROPConnector import ThermodynamicPoint as TP
import numpy as np
import matplotlib.pyplot as plt

# %%------------   DATA                        -----------------------------------------------------------> #
n = 100
P_in = 997086           # [Pa]
x_in = 0                # [-]
P_out = 736759          # [Pa]
h_out = 269200.98973    # [J/kgK]
m_dot_s = 0.050794957   # [kg/s]
H_s = 0.0005            # [m]
d_trot = 0.008

base_point = TP(["R1234ze"], [1], unit_system="MASS BASE SI")

input_point = base_point.duplicate()
input_point.set_variable("P", P_in)
input_point.set_variable("x", x_in)

H0 = input_point.get_variable("h")
S0 = input_point.get_variable("S")

output_point = base_point.duplicate()
output_point.set_variable("P", P_out)
output_point.set_variable("h", h_out)

H1 = output_point.get_variable("h")
S1 = output_point.get_variable("S")

# %%------------   CALCULATIONS                        -----------------------------------------------------------> #

Hh = np.linspace(H0, H1, n)
Ph = np.linspace(P_in, P_out, n)
Ss = np.linspace(S0, S0, n)
rhov_throat = 0.

i_max = 0.
v_sound = 0.
eta_nozzle = 0.85

v = np.zeros(n)
rhov = np.zeros(n)

m_dot_arr = np.zeros(n)

int_points = list()
for i in range(n):
    int_points.append(base_point.duplicate())

for i in range(n):

    int_points[i].set_variable("P", Ph[i])
    int_points[i].set_variable("s", Ss[i])

    h_curr = H0 - eta_nozzle * (H0 - int_points[i].get_variable("h"))

    int_points[i].set_variable("P", Ph[i])
    int_points[i].set_variable("h", h_curr)

    v[i] = np.sqrt(2 * (H0 - h_curr))
    rhov[i] = v[i] * int_points[i].get_variable("rho")

    if rhov[i] > rhov_throat:

        i_max = i
        v_sound = v[i]
        rhov_throat = rhov[i]

    A_throat = d_trot * H_s
    m_dot_s = A_throat * rhov_throat
    h_out = H0 - (v_sound ** 2) / 2

    m_dot_arr[i] = m_dot_s

output_point.set_variable("P", P_out)
output_point.set_variable("h", h_out)

# %%------------   PLOT                        -----------------------------------------------------------> #

fig = plt.subplots()
plt.plot(Ph, m_dot_arr)

plt.show()

