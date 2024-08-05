# %%-------------------------------------   IMPORT CLASSES   --------------------------------------------------------> #

import numpy as np
import matplotlib.pyplot as plt

# %%--------------------------------------   DATA COLLECT   ---------------------------------------------------------> #

disk_gap = [0.4, 0.6, 1, 1.25]

CFD_vr_plot = [0.2414, 0.2529, 0.2667, 0.2747]
CFD_vr_table = [0.21807, 0.22841, 0.23903, 0.24668]
EES_vr = [0.2667, 0.2724, 0.2747, 0.2805]
PY_vr = [0.2132, 0.2183, 0.2239, 0.2292]

CFD_torque = [0.00299, 0.00364, 0.00539, 0.00662]
EES_torque = [0.00403, 0.00591, 0.00864, 0.01104]
PY_torque = [0.004441, 0.006090, 0.008846, 0.011336]

# %%--------------------------------------   PLOT RESULTS   ---------------------------------------------------------> #

fig, axs = plt.subplots(2, 1, constrained_layout=True)

axs[0].plot(disk_gap, EES_vr, label='EES Model Plot', linestyle='', marker='s', mec='darkred', mfc='white')
axs[0].plot(disk_gap, CFD_vr_plot, label='CFD Model Plot', linestyle='', marker='o', mec='darkblue', mfc='white')
axs[0].plot(disk_gap, CFD_vr_table, label='CFD Model Table', linestyle='', marker='^', mec='darkblue', mfc='white')
axs[0].plot(disk_gap, PY_vr, label='PY Model', linestyle='', marker='x', mec='black', mfc='black')
axs[0].set_xlim([0.2, 1.3])
axs[0].set_ylim([0.2, 0.4])
axs[0].set_xlabel('Disk Gap [mm]')
axs[0].set_ylabel('Velocity Ratio [-]')
axs[0].set_title('Velocity Ratio')
axs[0].legend(loc='upper left')
axs[0].grid()

axs[1].plot(disk_gap, EES_torque, label='EES Model', linestyle='', marker='s', mec='darkred', mfc='white')
axs[1].plot(disk_gap, CFD_torque, label='CFD Model', linestyle='', marker='o', mec='darkblue', mfc='white')
axs[1].plot(disk_gap, PY_torque, label='PY Model', linestyle='', marker='x', mec='black', mfc='black')
axs[1].set_xlim([0.2, 1.6])
axs[1].set_ylim([0, 0.012])
axs[1].set_xlabel('Disk Gap [mm]')
axs[1].set_ylabel('Torque [N*m]')
axs[1].set_title('Torque')
axs[1].legend()
axs[1].grid()

plt.show()
