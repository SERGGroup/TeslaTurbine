# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.constants import CALCULATION_FOLDER, os
import matplotlib.pyplot as plt
import numpy as np

# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
result_folder = os.path.join(

    CALCULATION_FOLDER, "1 - Federico Schneider Design", "res"

)

filename = os.path.join(result_folder, "opt_results.npy")
results = np.load(filename)

filename = os.path.join(result_folder, "outlet_pressure_mesh.npy")
p_out = np.load(filename)

filename = os.path.join(result_folder, "throat_widths.npy")
tw = np.load(filename)

filename = os.path.join(result_folder, "turb_diameters_mesh.npy")
Diameters = np.load(filename)

res_shape = results.shape
n = res_shape[-1]
n_tw = res_shape[1]

# %%-------------------------------------   PLOT LINES                          -------------------------------------> #
fig, axs = plt.subplots(1, 3, dpi=300)
fig.set_size_inches(10, 5)
cmap = plt.get_cmap('viridis')
norm = plt.Normalize(0, 1)

i = 1
for j in range(n):

    color = cmap(norm(j / (n - 1)))

    for k in range(len(axs)):

        axs[k].plot(tw[:, j], results[k, i, :, j], color=color, zorder=1)
        axs[k].set_xscale("log")

        if k < len(axs) - 1:
            axs[k].set_yscale("log")

plt.tight_layout(pad=2)
plt.subplots_adjust(hspace=0)
plt.show()