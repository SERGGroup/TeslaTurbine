# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.single_phase import SPRotor, SPTeslaGeometry, SPTeslaOptions, SPStatorMil, SPStator
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from tqdm import tqdm

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

# Input Constraints Data
P_in = 10100000  # [Pa]
P_out = 3900000         # [Pa]
T_in_c = 94  # [°C]
T_in = T_in_c + 273.15  # [K]

m_rif = 0.157         # [kg/s]

setpoints = 50
rpm_array = np.linspace(5000, 30000, setpoints)

output_array = np.empty((len(rpm_array), 12))

# %%------------   CALCULATION                             ----------------------------------------------------------> #

for i in tqdm(range(setpoints)):

    # INITIALIZING TURBINE CONDITIONS
    curr_geometry = SPTeslaGeometry()
    curr_options = SPTeslaOptions()
    curr_options.rotor.integr_variable = 0.04  # [-]
    curr_options.stator.metastability_check = False
    curr_options.rotor.sp_check = True
    curr_options.rotor.n_rotor = 4000

    # Main design Parameters
    curr_geometry.rotor.b_channel = 0.00005
    curr_geometry.rotor.d_ratio = 3  # [m]
    curr_geometry.d_main = 0.15  # [m]

    tt = BaseTeslaTurbine("CarbonDioxide", curr_geometry, curr_options, stator=SPStatorMil, rotor=SPRotor)

    tt.rotor.gap_losses_control = False
    tt.rotor.rpm = rpm_array[i]

    Z_stat = 1
    tt.geometry.stator.Z_stat = Z_stat

    throat_width = 0.00115
    tt.geometry.throat_width = throat_width

    tt.geometry.disc_thickness = 0.0001
    tt.geometry.rotor.n_discs = 1
    tt.geometry.H_s = 0.00015
    tt.m_dot_tot = m_rif
    tt.geometry.n_channels = tt.n_packs

    alpha_in = 89  # [°]
    tt.geometry.alpha1 = alpha_in

    tt.points[0].set_variable("P", P_in)
    tt.points[0].set_variable("T", T_in)
    tt.stator.stator_eff = 0.9
    tt.rotor.gap_losses_control = True

    tt.P_in = P_in
    tt.P_out = P_out
    tt.T_in = T_in

    tt.iterate_pressure()
    tt.evaluate_performances()
    rotor_array = tt.rotor.get_rotor_array()

    output_array[i, 0] = tt.rotor.rpm
    output_array[i, 1] = tt.Eta_tesla_tt2
    output_array[i, 2] = tt.work
    output_array[i, 3] = tt.power
    output_array[i, 4] = tt.work2
    output_array[i, 5] = tt.stator.m_dot_s
    output_array[i, 6] = tt.n_packs
    output_array[i, 7] = tt.rotor.first_speed.vt / (tt.rotor.omega * tt.geometry.rotor.r_out)
    output_array[i, 8] = tt.stator.out_speed
    output_array[i, 9] = (tt.static_points[2].get_variable("h") - tt.static_points[3].get_variable("h")) / (
                tt.static_points[0].get_variable("h") - tt.static_points[3].get_variable("h"))
    output_array[i, 10] = tt.Eta_tesla_ss
    output_array[i, 11] = tt.Eta_tesla_ts


# %%------------   PLOT RESULTS                            ----------------------------------------------------------> #

fig, axs = plt.subplots(1,2, figsize=(8.6, 4))

line1, = axs[0].plot(output_array[:, 0], output_array[:, 3], label='Power per Channel', color='Darkblue', linewidth=1.3)
axs[0].set_xlabel('Rotational Speed [rpm]')
axs[0].set_ylabel('Power [W]')

x_start, x_end = 5000, 30000
y1_start, y1_end = 40, 160

axs[0].set_xlim(x_start, x_end)
axs[0].set_ylim(y1_start, y1_end)

rect1 = patches.Rectangle(
    (x_start, y1_start), 10000, y1_end - y1_start, linewidth=1, edgecolor="black",
    facecolor="orange", alpha=0.3
)

axs[0].add_patch(rect1)

ax2 = axs[0].twinx()
line2, = ax2.plot(output_array[:, 0], output_array[:, 7], ls='--', label='Tip Ratio', color='black', linewidth=1.3, zorder=1)
ax2.set_ylabel('Tip Ratio [-]')
ax2.set_ylim(0, 6)

handles, labels = axs[0].get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()
combined_handles = handles + handles2
combined_labels = labels + labels2
axs[0].legend(handles=combined_handles, labels=combined_labels, loc='upper right', edgecolor='black', facecolor='white')
axs[0].grid()
axs[0].set_title('a)', fontsize=10)

axs[0].annotate('Generator Operating Area', xy = (5000, 0.8), fontsize = 9, xytext = (16000, 0.2), bbox = dict(fc = 'white', boxstyle="round,pad=0.3", ec = 'black'), arrowprops = dict(fc = 'black', width = 0.05, headwidth = 4, headlength = 4), zorder = 1)

line3, = axs[1].plot(output_array[:, 0], output_array[:, 11], label='Efficiency', color='Darkred', linewidth=1.3)
axs[1].set_xlabel('Rotational Speed [rpm]')
axs[1].set_ylabel('Total-to-Static Efficiency [-]')

y2_start, y2_end = 0.2, 0.8

axs[1].set_ylim(y2_start, y2_end)
axs[1].set_xlim(5000, 30000)

rect2 = patches.Rectangle(
    (x_start, y2_start), 10000, y2_end - y2_start, linewidth=1, edgecolor="orange",
    facecolor="orange", alpha=0.3
)

axs[1].add_patch(rect2)

ax1 = axs[1].twinx()
line4, = ax1.plot(output_array[:, 0], output_array[:, 7], ls='--', label='Tip Ratio', color='Black', linewidth=1.3)
ax1.set_ylabel('Tip Ratio [-]')
ax1.set_ylim(0,6)

handles1, labels1 = axs[1].get_legend_handles_labels()
handles3, labels3 = ax1.get_legend_handles_labels()
combined_handles1 = handles1 + handles3
combined_labels1 = labels1 + labels3
axs[1].legend(handles=combined_handles1, labels=combined_labels1, loc='upper right', edgecolor='black', facecolor='white')
axs[1].annotate('Generator Operating Area', xy = (12000, 0.28), fontsize = 9, xytext = (16000, 0.22),
             bbox = dict(fc = 'white', boxstyle="round,pad=0.3", ec = 'black'),
             arrowprops = dict(fc = 'black', width = 0.05, headwidth = 4, headlength = 4)
             )
axs[1].grid()
axs[1].set_title('b)', fontsize=10)

plt.tight_layout()
axs[0].annotate('Generator Operating Area', xy = (12000, 56), fontsize = 9, xytext = (16000, 44),
             bbox = dict(fc = 'white', boxstyle="round,pad=0.3", ec = 'black'),
             arrowprops = dict(fc = 'black', width = 0.05, headwidth = 4, headlength = 4)
             )
plt.show()

