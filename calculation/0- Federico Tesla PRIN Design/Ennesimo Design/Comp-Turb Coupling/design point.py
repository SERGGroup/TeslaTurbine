# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.single_phase import SPRotor, SPTeslaGeometry, SPTeslaOptions, SPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

# Input Constraints Data
P_out = 3900000         # [Pa]
T_sat_c = 4.312 # [°C]
T_sat = T_sat_c + 273.15  # [K]
T_in_c = 94  # [°C]
T_in = T_in_c + 273.15  # [K]

m_rif = 0.1584             # [kg/s]

perc = np.linspace(0.7,1, 50)


# %%------------   CALCULATION                             ----------------------------------------------------------> #

# INITIALIZING OUTPUT VECTORS
inlet_pressure = list()
flowrate_ratio = list()
power = list()
efficiency = list()

for j in tqdm(range(len(perc))):

    m_comp = 0.
    P_in = 0.
    n_it = 0

    # INITIALIZING ITERATIVE PROCEDURE
    P_up = 13000000
    P_down = 8000000

    n_max = 100

    for i in tqdm(range(n_max)):

        P_in = (P_up + P_down) / 2

        P_in_Bar = P_in / 100000
        # COMPRESSOR CURVE
        m_comp = 0.250850086344 + 0.00689092778088 * T_sat_c - 0.00095376084804 * P_in_Bar + 0.0000645108060168 * T_sat_c ** 2 - 0.00000537648636096 * P_in_Bar * T_sat_c + 0.00000281599522776 * P_in_Bar ** 2

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
        tt.rotor.rpm = 15000

        Z_stat = 1
        tt.geometry.stator.Z_stat = Z_stat

        throat_width = 0.00115
        tt.geometry.throat_width = throat_width

        tt.geometry.disc_thickness = 0.0001
        tt.geometry.rotor.n_discs = 1
        tt.geometry.H_s = 0.00015
        # tt.m_dot_tot = m_rif
        tt.n_packs = 33
        tt.geometry.n_channels = tt.n_packs

        alpha_in = 89  # [°]
        tt.geometry.alpha1 = alpha_in
        tt.stator.stator_eff = 0.9
        tt.rotor.gap_losses_control = True
        tt.points[0].set_variable("P", P_in)
        tt.points[0].set_variable("T", T_in)

        tt.P_in = P_in
        tt.P_out = P_out
        tt.T_in = T_in

        tt.iterate_pressure()
        tt.evaluate_performances()
        rotor_array = tt.rotor.get_rotor_array()

        if abs(tt.stator.m_dot_s * tt.n_packs - perc[j] * m_comp) < 0.0001:

            break

        else:

            if (tt.stator.m_dot_s * tt.n_packs - perc[j] * m_comp) > 0:
                P_up = P_in
            else:
                P_down = P_in

        n_it += 1

    power.append(tt.power * tt.n_packs)
    efficiency.append(tt.Eta_tesla_ts)
    flowrate_ratio.append(tt.m_dot_tot / m_comp)
    inlet_pressure.append(P_in)

    print(' ')
    print('The {}-th percentage calculation has converged after {} iterations.'.format(j + 1, n_it))
    print(' ')

# %%------------                        PLOT                        -------------------------------------------------> #
plt.plot(perc, inlet_pressure)
plt.show()

# %%------------             TRAJECTORY PLOT                        -------------------------------------------------> #
fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={"projection": "polar"}, figsize=(10, 6))

ax1.plot(rotor_array[:, 24] * np.pi / 180, rotor_array[:, 1], color='Darkblue', linewidth='1.5')
ax1.set_rmax(0.075)
ax1.set_rticks([tt.geometry.d_main / tt.geometry.rotor.d_ratio /2,
               (tt.geometry.d_main - tt.geometry.d_main / tt.geometry.rotor.d_ratio) / 2])
ax1.grid(True)

ax1.set_title("Absolute Reference Frame")

ax2.plot(rotor_array[:, 25] * np.pi / 180, rotor_array[:, 1], color='Darkred', linewidth='1.5')
ax2.set_rmax(0.075)
ax2.set_rticks([tt.geometry.d_main / tt.geometry.rotor.d_ratio /2,
               (tt.geometry.d_main - tt.geometry.d_main / tt.geometry.rotor.d_ratio) / 2])
ax2.grid(True)

ax2.set_title("Relative Reference Frame")

plt.show()