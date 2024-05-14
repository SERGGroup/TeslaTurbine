import concurrent.futures

from main_code.constants import CALCULATION_FOLDER, os
from main_code.sub_classes.multi_phase import TPRotor, TPTeslaGeometry, TPTeslaOptions, TPStatorMil
from main_code.base_classes import BaseTeslaTurbine
from scipy.optimize import minimize
import numpy as np
import math
import time
import gc

result_folder = os.path.join(

    CALCULATION_FOLDER, "1 - Federico Schneider Design", "res"

)

support_folder = os.path.join(result_folder, "support")

curr_geometry = TPTeslaGeometry()
curr_geometry.rotor.n_channels = 1  # [-]
curr_geometry.rotor.b_channel = 0.0029
curr_geometry.rotor.roughness = 0.0000005

curr_options = TPTeslaOptions()
curr_options.rotor.profile_rotor = True
curr_options.rotor.sp_check = False
curr_options.rotor.tp_epsilon_model = "sarti"

fluid = "R1234ze"
tesla_turbine = BaseTeslaTurbine("R1234ze", curr_geometry, curr_options, stator=TPStatorMil, rotor=TPRotor)

P_in = 997233         # [Pa]
T_in = 323.15         # [K]

tesla_turbine.points[0].set_variable("P", P_in)
tesla_turbine.points[0].set_variable("T", T_in)
tesla_turbine.stator.stator_eff = 0.81
tesla_turbine.rotor.omega = 177.9
tesla_turbine.rotor.gap_losses_control = True

P_ref = 427308      # [Pa]
d_ref = 0.5         # [m]
tw_ref = 0.003      # [m]
n = 3
tw_s = np.linspace(0.5, 1.5, 3) * tw_ref
n_tw = len(tw_s)

P_out_s_list = np.linspace(0.98, 1.02, n) * P_ref
D_turb_list = np.linspace(0.9, 1.1, n) * d_ref
P_out_s, D_turb = np.meshgrid(P_out_s_list, D_turb_list)

base_filename = "res_{i}_{j}.npy"
reporting_points = [0.05, 0.25, 0.5, 0.75, 1]

n_opt_params = 3

def turbine_calculation(i, j, k):

    tesla_turbine.P_in = P_in
    tesla_turbine.P_out = P_out_s[j, k]
    curr_geometry.stator.d_int = D_turb[j, k]
    curr_geometry.throat_width = tw_s[i]

    tesla_turbine.solve_with_stator_outlet_pressure(P_out_s[j, k])
    rotor_array = tesla_turbine.rotor.get_rotor_array()
    tesla_turbine.evaluate_performances()


def evaluate(i, j):

    rep_pos = 0
    curr_value = i * n + j + 1
    max_digits = math.ceil(np.log10(n * n_tw))
    current_name = "{{:{}d}} of {{}}".format(max_digits).format(curr_value, n * n_tw)

    print("{} - Started".format(current_name))

    sub_start = time.time()

    sub_results = np.empty([n_opt_params, n])
    sub_results[:] = np.nan

    for k in range(n):
        def opt_func(x, opt_value=(0,)):

            dv_perc = x[0]
            tesla_turbine.rotor.dv_perc = dv_perc
            turbine_calculation(i=i, j=j, k=k)
            print("{}-{}".format(dv_perc, tesla_turbine.eta_tt))
            if opt_value[0] == 0:

                return -tesla_turbine.eta_tt

            elif opt_value[0] == 1:

                return -tesla_turbine.work

            else:

                return -tesla_turbine.power

        x0 = [1, 1, 1]

        for q in range(3):

            res = minimize(opt_func, np.array((x0[q])), args=[q])
            opt_res = -res.fun

            sub_results[q, k] = opt_res

        gc.collect()

        if k/n > reporting_points[rep_pos]:

            time_elapsed = time.time() - sub_start
            time_remaining = time_elapsed * (1 - reporting_points[rep_pos]) / reporting_points[rep_pos]

            print("{} - {:.0f}% done! -> {:.0f}s ({:.0f}min) elapsed, {:.0f}s ({:.0f}min) to-go".format(

                current_name, reporting_points[rep_pos]*100,
                time_elapsed, time_elapsed/60,
                time_remaining, time_remaining/60

            ))

            rep_pos += 1

    filename = os.path.join(support_folder, base_filename.format(i=i, j=j))
    np.save(filename, sub_results)

    time_elapsed = time.time() - sub_start

    print("{} - Done! -> {:.0f}s ({:.0f}min) elapsed".format(current_name, time_elapsed, time_elapsed / 60))

    return time_elapsed


if __name__ == '__main__':

    evaluate(0, 0)

    # print("Calculation Started!")
    #
    # with concurrent.futures.ProcessPoolExecutor() as executor:
    #
    #     future_list = list()
    #
    #     for i in range(n_tw):
    #
    #         for j in range(n):
    #
    #             time.sleep(0.5)
    #             future_list.append(executor.submit(evaluate, i, j))
    #
    #     elapsed_time_list = list()
    #
    #     for f in concurrent.futures.as_completed(future_list):
    #
    #         try:
    #             elapsed_time_list.append(f.result())
    #
    #         except:
    #             elapsed_time_list.append(np.nan)
    #
    # print("Calculation Completed!!")
    # print("Saving Results ...")
    # time.sleep(2)
    #
    # res_shape = np.append([n_opt_params, n], D_turb.shape)
    # results = np.empty(res_shape)
    # results[:] = np.nan
    #
    # for i in range(n_tw):
    #
    #     for j in range(n):
    #
    #         filename = os.path.join(support_folder, base_filename.format(i=i, j=j))
    #
    #         if os.path.exists(filename):
    #
    #             results[:, i, j, :] = np.load(filename)
    #             # os.remove(filename)
    #
    # filename = os.path.join(result_folder, "opt_results.npy")
    # np.save(filename, results)
    #
    # filename = os.path.join(result_folder, "outlet_pressure_mesh.npy")
    # np.save(filename, P_out_s)
    #
    # filename = os.path.join(result_folder, "turb_diameters_mesh.npy")
    # np.save(filename, D_turb)
    #
    # filename = os.path.join(result_folder, "throat_widths.npy")
    # np.save(filename, tw_s)
    #
    # print("Calculation Completed!!!")