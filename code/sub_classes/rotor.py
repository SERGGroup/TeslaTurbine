from code.tesla_turbine_class import TeslaTurbine
import numpy as np

class Rotor:

    def __init__(self, main_turbine: TeslaTurbine):

        self.main_turbine = main_turbine
        self.options = self.main_turbine.options

        self.input_point = self.main_turbine.points[2]
        self.output_point = self.main_turbine.points[3]

        self.rotor_points = list()

    def solve(self):

        first_step = RotorStep(self, r_max)
        new_step = first_step

        for i in range(self.options.n_rotor):

            if self.options.profile_rotor:
                self.rotor_points.append(new_step)

            new_step = new_step.update()
    def get_rotor_array(self):

        rotor_array = np.empty((self.options.n_rotor, 30))
        rotor_array[:] = np.nan

        if self.options.profile_rotor:

            for i in range(len(self.rotor_points)):

                curr_rotor = self.rotor_points[i]

                rotor_array[i, 0] = i

                rotor_array[i, 1] = curr_rotor.r
                rotor_array[i, 2] = curr_rotor.u

                rotor_array[i, 3] = curr_rotor.vt
                rotor_array[i, 4] = curr_rotor.vr
                rotor_array[i, 5] = curr_rotor.wt
                rotor_array[i, 6] = curr_rotor.wr

                rotor_array[i, 7] = curr_rotor.alpha
                rotor_array[i, 8] = curr_rotor.beta

                rotor_array[i, 9] = curr_rotor.v
                rotor_array[i, 10] = curr_rotor.w

                rotor_array[i, 11] = curr_rotor.thermo_point.get_variable("P")
                rotor_array[i, 12] = curr_rotor.thermo_point.get_variable("h")
                rotor_array[i, 13] = curr_rotor.thermo_point.get_variable("T")
                rotor_array[i, 14] = curr_rotor.thermo_point.get_variable("s")
                rotor_array[i, 15] = curr_rotor.thermo_point.get_variable("rho")
                rotor_array[i, 16] = curr_rotor.thermo_point.get_variable("mu")

                rotor_array[i, 17] = curr_rotor.Ma_R_a
                rotor_array[i, 18] = curr_rotor.Ma_R_r

                rotor_array[i, 19] = curr_rotor.Re_a
                rotor_array[i, 20] = curr_rotor.thermo_point.get_variable("quality")
                rotor_array[i, 21] = curr_rotor.cosy

                rotor_array[i, 24] = curr_rotor.theta
                rotor_array[i, 25] = curr_rotor.gamma
                rotor_array[i, 26] = curr_rotor.admr
                rotor_array[i, 27] = curr_rotor.theta(90)
                rotor_array[i, 28] = curr_rotor.theta(180)
                rotor_array[i, 29] = curr_rotor.theta(270)

        return rotor_array


class RotorStep:

    def __init__(self, main_rotor: Rotor, r):

        self.r = r
        self.main_rotor = main_rotor
        self.thermo_point = self.main_rotor.input_point.duplicate()
