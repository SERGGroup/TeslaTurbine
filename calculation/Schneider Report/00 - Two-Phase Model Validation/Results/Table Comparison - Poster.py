# %%-------------------------------------   IMPORT CLASSES   --------------------------------------------------------> #
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# %%--------------------------------------   DATA COLLECT   ---------------------------------------------------------> #

disk_gap = [0.8, 1.2, 2, 2.5]

CFD_vr_table = [0.21807, 0.22841, 0.23903, 0.24668]
PY_vr = [0.2132, 0.2183, 0.2239, 0.2292]

vr_error = [-2.23, -4.43, -6.33, -7.09]

CFD_pr = [1.086685, 1.088052, 1.083509, 1.087670]
PY_pr = [1.118288, 1.120646, 1.123126, 1.124680]

pr_error = [2.91, 3.00, 3.66, 3.40]

# %%--------------------------------------   PLOT RESULTS   ---------------------------------------------------------> #

fig, axs = plt.subplots(2, 2, constrained_layout=True, figsize=[8, 4])

def format_func(value, tick_number):
    return f"{value:.2g}"

axs[0, 0].plot(disk_gap, CFD_vr_table, label='CFD Model', linestyle='', marker='D', mec='black', mfc='white')
axs[0, 0].plot(disk_gap, PY_vr, label='Python Model', linestyle='', marker='^', mec='navy', mfc='navy')
axs[0, 0].set_xlabel('Disk Gap / Reference Disk Gap [-]')
axs[0, 0].set_ylabel('Velocity Ratio [-]')
axs[0, 0].set_xlim(0.4, 3.2)
axs[0, 0].set_ylim(0.2, 0.3)
axs[0, 0].legend(loc='upper left', edgecolor='black', facecolor= 'white')
axs[0, 0].grid()

axs[1, 0].bar(disk_gap, vr_error, width=0.25, color='navy', zorder=3)
axs[1, 0].set_xlabel('Disk Gap / Reference Disk Gap [-]')
axs[1, 0].set_ylabel('Relative Error [%]')
axs[1, 0].set_xlim(0.4, 3.2)
axs[1, 0].set_ylim(-8, 8)
axs[1, 0].grid(axis='y', zorder=0)

axs[1, 0].yaxis.set_major_formatter(ticker.FuncFormatter(format_func))

axs[0, 1].plot(disk_gap, CFD_pr, label='CFD Model', linestyle='', marker='D', mec='black', mfc='white')
axs[0, 1].plot(disk_gap, PY_pr, label='Python Model', linestyle='', marker='^', mec='navy', mfc='navy')
axs[0, 1].set_xlabel('Disk Gap / Reference Disk Gap [-]')
axs[0, 1].set_ylabel('Pressure Ratio [-]')
axs[0, 1].set_xlim(0.4, 3.2)
axs[0, 1].set_ylim(1, 1.4)
axs[0, 1].legend(loc='upper left', edgecolor='black', facecolor= 'white')
axs[0, 1].grid()

axs[1, 1].bar(disk_gap, pr_error, width=0.25, color='navy',zorder=3)
axs[1, 1].set_xlabel('Disk Gap / Reference Disk Gap [-]')
axs[1, 1].set_ylabel('Relative Error [%]')
axs[1, 1].set_xlim(0.4, 3.2)
axs[1, 1].set_ylim(-8, 8)
axs[1, 1].grid(axis='y', zorder=0)

axs[1, 1].yaxis.set_major_formatter(ticker.FuncFormatter(format_func))

plt.show()