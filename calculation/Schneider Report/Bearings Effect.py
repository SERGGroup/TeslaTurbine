# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.sub_classes.multi_phase import TPRotor, TPTeslaGeometry, TPTeslaOptions, TPStatorMil
from main_code.base_classes import BaseTeslaTurbine
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from main_code.base_classes.support.bearing_losses import BearingLoss

# %%------------   MAIN INPUT DATA                         ----------------------------------------------------------> #

n_setpoints = 40
dv_perc_list = np.linspace(-0.8, 0.2, n_setpoints)
TW_ref = 0.003  # [m]
TW_list = np.linspace(1, 10, n_setpoints) * TW_ref
dv_perc, TW = np.meshgrid(dv_perc_list, TW_list)

curr_options = TPTeslaOptions()
curr_options.rotor.profile_rotor = True
curr_options.rotor.sp_check = False
curr_options.rotor.tp_epsilon_model = "sarti"
# curr_options.stator.metastability_check = True

P_in = 997233  # [Pa]
x_in = 0  # [-]
P_out = 427308  # [Pa]
m_refr = 3.336  # [kg/s]

output_array_list = list()
rotor_array_max_efficiency_list = list()
rotor_array_max_power_list = list()

Eta = np.empty(len(dv_perc_list) * len(TW_list))
Power = np.empty(len(dv_perc_list) * len(TW_list))

Eta_arr = np.empty(len(dv_perc_list) * len(TW_list))
Power_arr = np.empty(len(dv_perc_list) * len(TW_list))
RPM_arr = np.empty(len(dv_perc_list) * len(TW_list))

output_array = np.empty((len(range(n_setpoints)) * len(range(n_setpoints)), 13))
# %%------------   CALCULATION                             ----------------------------------------------------------> #

for i in tqdm(range(n_setpoints)):

    print(' ')
    print(' ')
    print('Starting of {}-th dv_perc Calculation'.format(i + 1))

    j_0 = i * n_setpoints

    for j in tqdm(range(n_setpoints)):

        # Geometric Parameters
        curr_geometry = TPTeslaGeometry()
        curr_geometry.d_main = 0.4
        curr_geometry.rotor.b_channel = 0.0005  # [m]
        curr_geometry.disc_thickness = 0.0008  # [m]
        curr_geometry.alpha1 = 85  # [°]
        curr_geometry.stator.Z_stat = 4  # [-]

        curr_geometry.rotor.gap = 0.001             # [m]

        curr_geometry.throat_width = TW[i, j]
        curr_geometry.rotor.n_discs = 1             # [-]
        curr_geometry.rotor.d_ratio = 3.5           # [-]
        curr_geometry.rotor.n_packs = 1             # [-]
        curr_geometry.rotor.roughness = 0.0000005   # [m]

        single_hs = ( curr_geometry.rotor.n_discs - 1) * curr_geometry.disc_thickness + curr_geometry.rotor.n_discs * curr_geometry.rotor.b_channel
        curr_geometry.H_s = single_hs * curr_geometry.rotor.n_packs

        # Solving Parameters
        curr_options = TPTeslaOptions()
        curr_options.rotor.integr_variable = 0.03
        curr_options.rotor.n_rotor = 400

        new_turbine = BaseTeslaTurbine('R1234ze', curr_geometry, curr_options, stator=TPStatorMil, rotor=TPRotor)
        new_turbine.m_dot_tot = m_refr

        new_turbine.points[0].set_variable("P", P_in)
        new_turbine.points[0].set_variable("x", x_in)
        new_turbine.stator.stator_eff = 0.9
        new_turbine.rotor.gap_losses_control = True

        new_turbine.options.stator.metastability_check = True
        new_turbine.options.rotor.profile_rotor = True
        new_turbine.options.rotor.sp_check = False
        new_turbine.options.rotor.tp_epsilon_model = "sarti"


        new_turbine.rotor.dv_perc = dv_perc[i,j]
        # new_turbine.rotor.rpm = 1400
        new_turbine.m_dot_tot = m_refr

        new_turbine.P_in = P_in
        new_turbine.P_out = P_out
        new_turbine.iterate_pressure()
        # new_turbine.evaluate_performances()
        rotor_array = new_turbine.rotor.get_rotor_array()

        new_BL = BearingLoss()

        # Geometric Parameters
        new_BL.volume = new_turbine.volume
        new_BL.d = (new_turbine.geometry.d_main / new_turbine.geometry.rotor.d_ratio)   # Internal Bearing Diameter
        new_BL.D = new_BL.d + 0.01                                                      # External Bearing Diameter
        new_BL.dm = (new_BL.d + new_BL.D) / 2

        # Operational Parameters
        new_BL.nu = 100                                                                 # Lubricating Oil Viscosity in cst (or mm2/s)
        new_BL.n = new_turbine.rotor.rpm

        # Thermodynamic Parameters
        new_BL.m_dot = new_turbine.rotor.m_dot_ch * new_turbine.geometry.rotor.n_discs * new_turbine.n_packs
        new_BL.p = rotor_array[-1, 11]
        new_BL.v = rotor_array[-1, 9]
        new_BL.A_out = np.pi * (new_turbine.geometry.d_main / new_turbine.geometry.rotor.d_ratio) ** 2 / 4

        new_BL.evaluate_bearing_losses()

        new_turbine.evaluate_bearing_performances(new_BL.M_lost)

        Eta[j + j_0] = new_turbine.Eta_tesla_ss
        Power[j + j_0] = new_turbine.power * new_turbine.geometry.rotor.n_channels

        output_array[j + j_0, 0] = new_turbine.volume
        output_array[j + j_0, 1] = new_turbine.Eta_tesla_ss
        output_array[j + j_0, 2] = new_turbine.work
        output_array[j + j_0, 3] = new_turbine.power
        output_array[j + j_0, 4] = new_turbine.rotor.rpm
        output_array[j + j_0, 5] = new_turbine.stator.m_dot_s
        output_array[j + j_0, 6] = new_turbine.n_packs
        output_array[j + j_0, 7] = new_turbine.static_points[1].get_variable("p")
        output_array[j + j_0, 8] = new_turbine.stator.out_speed
        output_array[j + j_0, 9] = new_BL.M_lost
        output_array[j + j_0, 10] = new_turbine.rotor.omega
        output_array[j + j_0, 11] = new_BL.M_lost / (new_turbine.work * new_turbine.n_packs) * 100
        output_array[j + j_0, 12] = new_turbine.Eta_bearings

# %%------------   PLOT RESULTS                            ----------------------------------------------------------> #

res1 = output_array[:, 1].reshape(dv_perc.shape)
res2 = output_array[:, 12].reshape(dv_perc.shape)
res3 = ((output_array[:, 3] * output_array[:, 6]) / (output_array[:, 0] * 1000)).reshape(dv_perc.shape)
res4 = ((output_array[:, 3] * output_array[:, 6] - output_array[:, 9] * output_array[:, 10]) / (output_array[:, 0] * 1000)).reshape(dv_perc.shape)

fig, axs = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

ETA = axs[0, 0].contourf(dv_perc, TW, res1, levels=np.linspace(0,0.45, 15))
ETA_B = axs[0, 1].contourf(dv_perc, TW, res2, levels=np.linspace(0,0.45, 15))
SPEC_POWER = axs[1, 0].contourf(dv_perc, TW, res3, levels=np.linspace(0,350, 10))
SPEC_POWER_B = axs[1, 1].contourf(dv_perc, TW, res4, levels=np.linspace(0,350, 10))

axs[0, 0].set(xlabel='dv_perc [-]', ylabel='Throat Width [m]', title='Efficiency - No Bearings [-]')
axs[0, 1].set(xlabel='dv_perc [-]', ylabel='Throat Width [m]', title='Efficiency - Bearings [-]')
axs[1, 0].set(xlabel='dv_perc [-]', ylabel='Throat Width [m]', title='Specific Power - No Bearings [kW / m3]')
axs[1, 1].set(xlabel='dv_perc [-]', ylabel='Throat Width [m]', title='Specific Power - Bearings [kW / m3]')

CB1 = fig.colorbar(ETA, shrink=0.8, ax=axs[0, 0], aspect=30)
CB2 = fig.colorbar(ETA_B, shrink=0.8, ax=axs[0, 1], aspect=30)
CB3 = fig.colorbar(SPEC_POWER, shrink=0.8, ax=axs[1, 0], aspect=30)
CB4 = fig.colorbar(SPEC_POWER_B, shrink=0.8, ax=axs[1, 1], aspect=30)

plt.suptitle("Performance at D = 0.4 m")

plt.show()

# %%------------   PLOT RESULTS                            ----------------------------------------------------------> #

res1 = (100 * (output_array[:, 1] - output_array[:, 12])).reshape(dv_perc.shape)
res2 = ((100 * output_array[:, 9] * output_array[:, 10]) / (output_array[:, 3] * output_array[:, 6])).reshape(dv_perc.shape)

fig1, axs1 = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)

ETA = axs1[0].contourf(dv_perc, TW, res1)
POWER = axs1[1].contourf(dv_perc, TW, res2)

axs1[0].set(xlabel='dv_perc [-]', ylabel='Throat Width [m]', title='Efficiency Loss [%]')
axs1[1].set(xlabel='dv_perc [-]', ylabel='Throat Width [m]', title='Power Loss / Turbine Power [%]')

CB5 = fig1.colorbar(ETA, shrink=0.8, ax=axs1[0], aspect=30)
CB6 = fig1.colorbar(POWER, shrink=0.8, ax=axs1[1], aspect=30)

plt.suptitle("Bearings Effect")

plt.show()

