# %%-------------------------------------   IMPORT CLASSES   --------------------------------------------------------> #
import numpy as np
import matplotlib.pyplot as plt

# %%--------------------------------------   DATA COLLECT   ---------------------------------------------------------> #

disk_gap = [0.4, 0.6, 1, 1.25]

CFD_DeltaT = [3.52, 3.45, 3.55, 3.42]
PY_DeltaT = [4.823666, 4.914017, 5.003756, 5.066676]

CFD_vr_table = [0.21807, 0.22841, 0.23903, 0.24668]
PY_vr = [0.2132, 0.2183, 0.2239, 0.2292]

CFD_pr = [1.086685, 1.088052, 1.083509, 1.087670]
PY_pr = [1.118288, 1.120646, 1.123126, 1.124680]

# %%--------------------------------------   PLOT RESULTS   ---------------------------------------------------------> #

fig, axs = plt.subplots(3, 1, constrained_layout=True)

axs[0].plot(disk_gap, CFD_vr_table, label='CFD Model', linestyle='', marker='o', mec='darkblue', mfc='darkblue')
axs[0].plot(disk_gap, PY_vr, label='Python Model', linestyle='', marker='x', mec='darkred', mfc='darkred')
axs[0].set_xlabel('Disk Gap [mm]')
axs[0].set_ylabel('Velocity Ratio [-]')
axs[0].set_xlim(0.2, 1.6)
axs[0].set_ylim(0.2, 0.3)
axs[0].legend(loc='upper left', edgecolor='black', facecolor= 'white')
axs[0].grid()

axs[1].plot(disk_gap, CFD_DeltaT, label='CFD Model', linestyle='', marker='o', mec='darkblue', mfc='darkblue')
axs[1].plot(disk_gap, PY_DeltaT, label='Python Model', linestyle='', marker='x', mec='darkred', mfc='darkred')
axs[1].set_xlabel('Disk Gap [mm]')
axs[1].set_ylabel('Delta Temp [°C]')
axs[1].set_xlim(0.2, 1.6)
axs[1].set_ylim(2, 9)
axs[1].legend(loc='upper left', edgecolor='black', facecolor= 'white')
axs[1].grid()

axs[2].plot(disk_gap, CFD_pr, label='CFD Model', linestyle='', marker='o', mec='darkblue', mfc='darkblue')
axs[2].plot(disk_gap, PY_pr, label='Python Model', linestyle='', marker='x', mec='darkred', mfc='darkred')
axs[2].set_xlabel('Disk Gap [mm]')
axs[2].set_ylabel('Pressure Ratio [-]')
axs[2].set_xlim(0.2, 1.6)
axs[2].set_ylim(1, 1.4)
axs[2].legend(loc='upper left', edgecolor='black', facecolor= 'white')
axs[2].grid()

plt.show()