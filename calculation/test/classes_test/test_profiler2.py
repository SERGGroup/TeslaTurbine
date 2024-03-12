# %%------------   IMPORT CLASSES
import numpy as np
import matplotlib.pyplot as plt
from main_code.base_classes.support.geometry import BaseTeslaGeometry


# %%------------   DATA

n = 100
r_0 = 125
r_1 = 100
ax = np.zeros(n)

rad_axis = np.linspace(r_0, r_1, n)

geom = BaseTeslaGeometry()
x_ss, y_ss, y_ps, x_ps, m, alpha_out = geom.stator.get_coordinates(n, r_0, r_1, ax)

theta = np.linspace(0, 18, 100)
x_min = list()
y_min = list()
x_max = list()
y_max = list()

for theta in theta:
    x_min.append(r_0 * np.sin(theta * np.pi / 180))
    y_min.append(r_0 * np.cos(theta * np.pi / 180))
    x_max.append(r_1 * np.sin(theta * np.pi / 180))
    y_max.append(r_1 * np.cos(theta * np.pi / 180))

A_th, A_eff, throat_width = geom.stator.geometric_parameters(x_ss, y_ss, y_ps, x_ps, m, n)


# %%------------   PLOT

fig, axis = plt.subplots(1,2)

# Plotting the Geometric Profile of the Stator
axis[0].plot(x_ss, y_ss, label="Suction Side")
axis[0].plot(x_ps, y_ps, label="Pressure Side")
axis[0].plot(x_min, y_min, label="Internal Diameter")
axis[0].plot(x_max, y_max, label="External Diameter")

# Plotting the Evolution of the Sectional Area
axis[1].plot(rad_axis, A_th, label="Theoretical Area")
axis[1].plot(rad_axis, A_eff, label="Effective Area")

axis[0].set_xlabel("x [mm]")
axis[0].set_ylabel("y [mm]")

axis[1].set_xlabel("Rad. Direction [mm]")
axis[1].set_ylabel("Area [mm2]")

axis[0].set_ylim([80, 130])
axis[0].set_xlim([0, 30])
axis[1].set_xlim([100,125])

plt.show()




