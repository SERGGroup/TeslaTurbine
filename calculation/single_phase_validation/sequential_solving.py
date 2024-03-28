# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
import numpy as np
from main_code.sub_classes.single_phase import SPStator, SPRotor, SPTeslaGeometry, SPTeslaOptions
from main_code.base_classes import BaseTeslaTurbine
import matplotlib.pyplot as plt

# %%------------       INPUT DATA                         -----------------------------------------------------------> #
curr_geometry = SPTeslaGeometry()
curr_options = SPTeslaOptions()
curr_options.stator.iterate_phi = True
curr_geometry.stator.d_int = 0.2        # [m]
curr_options.rotor.integr_variable = 0.03
curr_options.rotor.n_rotor = 400

P_in = 22000000  # [Pa]
T_in = 423.15      # K
P_out = 21000000   # [Pa]

tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStator, rotor=SPRotor)
tt.points[0].set_variable("P", P_in)
tt.points[0].set_variable("T", T_in)

tt.static_points[1].set_variable("P", P_out)

tt.geometry.d_main = 0.2
tt.geometry.throat_width = 0.0004934
tt.geometry.H_s = 0.00094
tt.geometry.alpha1 = 85

tt.rotor.omega = 628.3185

tt.geometry.rotor.d_ratio = 2.5
tt.geometry.rotor.n_channels = 2

# %%------------     STATOR SOLVE                        -----------------------------------------------------------> #
tt.stator.solve()

# %%------------      ROTOR SOLVE                        -----------------------------------------------------------> #
tt.rotor.solve()
rotor_array = tt.rotor.get_rotor_array()
tt.evaluate_performances()
rpm = tt.rotor.rpm
P_out = rotor_array[-1, 11]
power = tt.power

# %%------------             TRAJECTORY PLOT                        -------------------------------------------------> #

fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

ax.plot(rotor_array[:, 24] *np.pi / 180, rotor_array[:, 1])
ax.set_rmax(0.1)
ax.set_rticks([0.025, 0.05, 0.075, 0.1])
ax.set_rlabel_position(-22.5)
ax.grid(True)

ax.set_title("theta_rel")

plt.show()

# %%------------             VELOCITY PLOT                        ---------------------------------------------------> #

plt.rcParams['figure.constrained_layout.use'] = True
fig1, (ax3, ax4) = plt.subplots(1, 2, figsize=(10, 6))

# Absolute Velocities
ax3.plot(rotor_array[:, 1], rotor_array[:, 3], label="Tangential Speed vt")
ax3.plot(rotor_array[:, 1], rotor_array[:, 4], label="Radial Speed vr")

ax3.set_xlabel("Radial Distance [m]", color = "Black", fontsize = 14)
ax3.set_ylabel("Velocity [m/s]", color = "Black", fontsize = 14)
ax3.set_xlim(0.04, 0.1)
ax3.set_ylim(0, 75)
ax3.legend(edgecolor = "Black", fontsize = 12)
ax3.grid()

# Relative Velocities
ax4.plot(rotor_array[:, 1], rotor_array[:, 5], label="Tangential Speed wt")
ax4.plot(rotor_array[:, 1], rotor_array[:, 6], label="Radial Speed wr")

ax4.set_xlabel("Radial Distance [m]", color = "Black", fontsize = 14)
ax4.set_ylabel("Velocity [m/s]", color = "Black", fontsize = 14)
ax4.set_xlim(0.04, 0.1)
ax4.set_ylim(0, 25)
ax4.legend(edgecolor = "Black", fontsize = 12)
ax4.grid()

plt.show()

# %%------------             VELOCITY PLOT 2                      ---------------------------------------------------> #

fig2 = plt.subplots()

plt.plot(rotor_array[:, 1], rotor_array[:, 2], label="Rotating Speed", c="Darkred")
plt.plot(rotor_array[:, 1], rotor_array[:, 9], label="Absolute Speed", c="Darkblue")
plt.plot(rotor_array[:, 1], rotor_array[:, 10], label="Relative Speed", c="Darkgreen")

plt.xlabel("Radial Distance [m]", color = "Black", fontsize = 14)
plt.ylabel("Velocity [m/s]", color = "Black", fontsize = 14)
plt.xlim(0.04, 0.1)
plt.ylim(0, 75)
plt.legend(edgecolor = "Black", fontsize = 12)
plt.grid()

plt.show()

# %%------------             THERMODYNAMIC PLOT                   ---------------------------------------------------> #

fig3 = plt.subplots()

plt.plot(rotor_array[:, 1], rotor_array[:, 11]/rotor_array[0, 11], label="Pressure", c="Darkred")
plt.plot(rotor_array[3:, 1], rotor_array[3:, 13]/rotor_array[3, 13], label="Temperature", c="Darkgreen")
plt.plot(rotor_array[3:, 1], rotor_array[3:, 15]/rotor_array[3, 15], label="Density", c="Orange")
plt.plot(rotor_array[3:, 1], rotor_array[3:, 12]/rotor_array[3, 12], label="Enthalpy", c="Darkblue")
plt.plot(rotor_array[3:, 1], rotor_array[3:, 14]/rotor_array[3, 14], label="Entropy", c="Olive")


plt.ylabel("Non-Dimensional Parameter [-]", fontsize = 14)
plt.xlabel("Radial Distance [m]", fontsize = 14)
plt.xlim(0.04, 0.1)
plt.ylim(0.93, 1.02)
plt.legend(edgecolor = "Black", fontsize = 12)
plt.grid()

plt.show()