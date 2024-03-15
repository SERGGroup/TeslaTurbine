# %%------------   IMPORT CLASSES
import numpy as np
import matplotlib.pyplot as plt
from main_code.base_classes.support.geometry import BaseTeslaGeometry, StatorGeometry

# %%------------   PROFILE COEFFICIENTS

r_0 = 125
r_int = 100
H_s = 0.6

# Parameters for Stator Profiling
t_u = 0.12  # [mm]
ch_u = 12.5  # [mm]
c_u = 0.19  # [mm]
b_u = 0.17  # [mm]
a_u = 0.1  # [mm]

t_b = 0.12  # [mm]
ch_b = 12.5  # [mm]
c_b = 0.28  # [mm]
b_b = 0.28  # [mm]
a_b = 0.1  # [mm]

r = np.linspace(r_0, r_int, 100)
theta = np.linspace(0, 15, 100)
ax = np.zeros(len(r))
up = np.zeros(len(r))
bottom = np.zeros(len(r))
midline = np.zeros(len(r))
m = np.zeros(len(r))
m_perp = np.zeros(len(r))
X_SS = np.zeros(len(r))
Y_SS = np.zeros(len(r))
X_PS = np.zeros(len(r))
Y_PS = np.zeros(len(r))
X_LM = np.zeros(len(r))
Y_LM = np.zeros(len(r))
A = np.zeros(len(r))
A_eff = np.zeros(len(r))

X_min = np.zeros(len(r))
Y_min = np.zeros(len(r))
X_max = np.zeros(len(r))
Y_max = np.zeros(len(r))

up_iniz = 12.5
bottom_iniz = 0.

for i in range(len(r)):
    ax[i] = r_0 - r[i]
    up[i] = up_iniz + 5 / 3 * t_u * ch_u * (a_u * (ax[i] / ch_u) ** 2 + b_u * (ax[i] / ch_u) ** 3 + c_u * (ax[i] / ch_u) ** 4)
    bottom[i] = bottom_iniz + 5 / 3 * t_b * ch_b * (a_b * (ax[i] / ch_b) ** 2 + b_b * (ax[i] / ch_b) ** 3 + c_b * (ax[i] / ch_b) ** 4)
    midline[i] = (up[i] + bottom[i]) / 2
    m[i] = 5 / 6 * t_u * (2 * a_u * (ax[i] / ch_u) + 3 * b_u * (ax[i] / ch_u) ** 2 + 4 * c_u * (ax[i] / ch_u) ** 3) + 5 / 6 * t_b * (2 * a_b * (ax[i] / ch_b) + 3 * b_b * (ax[i] / ch_b) ** 2 + 4 * c_b * (ax[i] / ch_b) ** 3)

    if m[i] == 0:
        m_perp[i] = -99999
    else:
        m_perp[i] = - 1 / m[i]

    X_SS[i] = bottom[i]
    Y_SS[i] = np.sqrt(r[i] ** 2 - bottom[i] ** 2)
    Y_PS[i] = np.sqrt(r[i] ** 2 - up[i] ** 2)
    X_PS[i] = np.sqrt(r[i] ** 2 - Y_PS[i] ** 2)
    X_LM[i] = (X_SS[i] + X_PS[i]) / 2
    Y_LM[i] = (Y_SS[i] + Y_PS[i]) / 2

    X_min[i] = r_0 * np.sin(theta[i] * np.pi / 180)
    Y_min[i] = r_0 * np.cos(theta[i] * np.pi / 180)
    X_max[i] = r_int * np.sin(theta[i] * np.pi / 180)
    Y_max[i] = r_int * np.cos(theta[i] * np.pi / 180)

    A[i] = np.sqrt((X_SS[i] - X_PS[i]) ** 2 + (Y_SS[i] - Y_PS[i]) ** 2) * H_s
    A_eff[i] = np.sqrt((X_SS[i] - X_PS[i]) ** 2 + (Y_SS[i] - Y_PS[i]) ** 2) * H_s * np.cos(np.arctan(m[i]))

Throat_width = np.min(A_eff) / H_s
alpha_out = (np.arcsin(X_LM[99] / r[99]) + np.arctan(m[99])) * 180 / np.pi

# %%------------   PLOT

# Stator Channel Shape
fig = plt.subplots()
plt.plot(X_SS, Y_SS)
plt.plot(X_PS, Y_PS)
plt.plot(X_min, Y_min)
plt.plot(X_max, Y_max)

# Throat Area
# fig2 = plt.subplots()
# plt.plot(r, A)
# plt.plot(r, A_eff)
plt.show()









