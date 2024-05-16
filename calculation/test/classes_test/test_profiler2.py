# %%------------   IMPORT CLASSES
import numpy as np
import matplotlib.pyplot as plt
from main_code.base_classes.support.geometry import BaseTeslaGeometry
from main_code.base_classes.support.options import StatorOptions


# %%------------   DATA

options = StatorOptions()
geom = BaseTeslaGeometry()

geom.d_main = 0.2

rad_axis = np.linspace(geom.d_main / 2 * geom.stator.d_ratio * 1000, geom.d_main / 2 * 1000, options.n_stator)
x_ss, y_ss, y_ps, x_ps, m, alpha, alpha_vert, dl = geom.stator.get_coordinates(options.n_stator)

chord = np.sum(dl)

theta = np.linspace(0, 18, 100)
x_min = list()
y_min = list()
x_max = list()
y_max = list()

for theta in theta:
    x_min.append(geom.d_main / 2 * 1000 * np.sin(theta * np.pi / 180))
    y_min.append(geom.d_main / 2 * 1000 * np.cos(theta * np.pi / 180))
    x_max.append(geom.d_main / 2 * 1000 * geom.stator.d_ratio * np.sin(theta * np.pi / 180))
    y_max.append(geom.d_main / 2 * 1000 * geom.stator.d_ratio * np.cos(theta * np.pi / 180))

A_th, A_eff = geom.stator.get_area(x_ss, y_ss, y_ps, x_ps, m, options.n_stator)

throat = geom.stator.throat_width
alpha_out = geom.stator.alpha_stat

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




