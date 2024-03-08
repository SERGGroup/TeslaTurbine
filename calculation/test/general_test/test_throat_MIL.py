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
H_s = 0.0005

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
Ss = np.linspace(S0, S0, n)
rhov_throat = 0.
v_sound = 0.

v = np.zeros(n)
rhov = np.zeros(n)

int_points = list()
for i in range(n):
    int_points.append(base_point.duplicate())

for i in range(n):
    int_points[i].set_variable("h", Hh[i])
    int_points[i].set_variable("s", Ss[i])
    v[i] = np.sqrt(2 * (H0 - Hh[i]))
    rhov[i] = v[i] * int_points[i].get_variable("rho")
    if rhov[i] > rhov_throat:
        v_sound = v[i]
        rhov_throat = rhov[i]

A_throat = m_dot_s / rhov_throat
D_throat = A_throat / H_s
# %%------------   PLOT                        -----------------------------------------------------------> #

fig = plt.subplots()

