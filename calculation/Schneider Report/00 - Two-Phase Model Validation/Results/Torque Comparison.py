# %%-------------------------------------   IMPORT CLASSES   --------------------------------------------------------> #

import numpy as np
import matplotlib.pyplot as plt

# %%--------------------------------------   DATA COLLECT   ---------------------------------------------------------> #

disk_gap = [0.4, 0.6, 1, 1.25]

CFD_torque = [0.00299, 0.00364, 0.00539, 0.00662]
EES_torque = [0.00403, 0.00591, 0.00864, 0.01104]
PY_torque = [0.004441, 0.006090, 0.008846, 0.011336]

# %%--------------------------------------   PLOT RESULTS   ---------------------------------------------------------> #

fig = plt.subplots(constrained_layout=True)

plt.plot(disk_gap, EES_torque, label='EES Model', linestyle='', marker='s', mec='black', mfc='black', mew=1.5, markersize = 4)
plt.plot(disk_gap, CFD_torque, label='CFD Model', linestyle=':', color='darkblue', marker='o', mec='darkblue',
         mfc='darkblue')
plt.plot(disk_gap, PY_torque, label='Python Model', linestyle='', marker='x', mec='darkred', mfc='darkred', mew=1.5)

plt.xlim(0.2, 1.6)
plt.ylim(0.002, 0.012)
plt.xlabel('Disk Gap [mm]', fontsize=12)
plt.ylabel('Torque [Nm]', fontsize=12)
plt.legend(loc='upper left', edgecolor='black', facecolor= 'white', fontsize=12)
plt.grid()

plt.show()
