from .support import Position, Speed
import numpy as np


class BaseRotorStep:

    def __init__(self, main_rotor, speed: Speed):

        self.main_rotor = main_rotor
        self.options = main_rotor.options

        self.speed = speed
        self.thermo_point = self.main_rotor.input_point.duplicate()

    @property
    def pos(self):

        return self.speed.pos

    @property
    def m_dot(self):

        return self.main_rotor.m_dot_ch

    def get_new_step(self, dr):

        # Evaluate main parameter variation
        dvt, dvr, dp, dh = self.get_variations(dr)

        # Update Position
        new_pos = self.speed.get_new_position(dr)

        # Update Speed
        new_speed = Speed(new_pos)
        vt_new = self.speed.vt + dvt
        vr_new = self.speed.vt + dvr
        new_speed.init_from_codes("vt", vt_new, "vr", vr_new)

        # Init a new step
        new_class = self.self_class()
        new_step = new_class(main_rotor=self.main_rotor, speed=new_speed)

        # Update Thermodynamic Point
        p_new = self.thermo_point.get_variable("P") + dp
        h_new = self.thermo_point.get_variable("H") + dh
        new_step.thermo_point.set_variable("P", p_new)
        new_step.thermo_point.set_variable("H", h_new)

        return new_step

    def get_variations(self, dr):

        return 0., 0., 0., 0.

    def __init_subclass__(cls, *args, **kwargs):

        def return_self_class(self):

            return type(self)

        cls.self_class = return_self_class


class BaseRotor:

    dv_perc = 0.

    def __init__(self, main_turbine, rotor_step: type(BaseRotorStep)):

        self.main_turbine = main_turbine
        self.geometry = self.main_turbine.geometry.rotor
        self.options = self.main_turbine.options.rotor

        self.input_point = self.main_turbine.points[2]
        self.output_point = self.main_turbine.points[3]

        self.__omega = 0.
        self.rotor_points = list()
        self.__rotor_step_cls = rotor_step

    def solve(self):

        self.__omega = self.main_turbine.stator.speed_out.vt / ((self.dv_perc + 1) * self.geometry.r_out)
        self.rotor_points = list()
        self.evaluate_gap_losses()

        first_pos = Position(self.geometry.r_out, self.__omega)
        first_speed = Speed(position=first_pos)
        first_speed.equal_absolute_speed_to(self.main_turbine.stator.speed_out)
        first_step = self.__rotor_step_cls(self, first_speed)

        new_step = first_step

        #
        # For detailed explanation on rotor discretization check:
        #
        #   "main_code/base_classes/other/rotor discretization explaination.xlsx"
        #

        dr_tot = self.geometry.dr_tot
        b = self.options.integr_variable * dr_tot
        a = np.power(dr_tot / b + 1, 1 / self.options.n_rotor)

        for i in range(self.options.n_rotor):

            if self.options.profile_rotor:
                self.rotor_points.append(new_step)

            dr = a ** i * (a - 1) * b
            new_step = new_step.get_new_step(dr)

        if self.options.profile_rotor:
            self.rotor_points.append(new_step)

    def evaluate_gap_losses(self):

        # No Losses Model
        self.main_turbine.points[1].copy_state_to(self.main_turbine.points[2])

    @property
    def m_dot_ch(self):

        return self.main_turbine.stator.m_dot_s / self.geometry.n_channels

    def get_rotor_array(self):

        rotor_array = np.empty((self.options.n_rotor, 30))
        rotor_array[:] = np.nan

        if self.options.profile_rotor:

            for i in range(len(self.rotor_points)):
                curr_rotor = self.rotor_points[i]

                rotor_array[i, 0] = i

                rotor_array[i, 1] = curr_rotor.pos.r
                rotor_array[i, 2] = curr_rotor.speed.u

                rotor_array[i, 3] = curr_rotor.speed.vt
                rotor_array[i, 4] = curr_rotor.speed.vr
                rotor_array[i, 5] = curr_rotor.speed.wt
                rotor_array[i, 6] = curr_rotor.speed.wr

                rotor_array[i, 7] = curr_rotor.speed.alpha
                rotor_array[i, 8] = curr_rotor.speed.beta

                rotor_array[i, 9] = curr_rotor.speed.v
                rotor_array[i, 10] = curr_rotor.speed.w

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

                rotor_array[i, 24] = curr_rotor.pos.theta_rel(0, "°")
                rotor_array[i, 25] = curr_rotor.pos.gamma_rel(0, "°")
                rotor_array[i, 26] = curr_rotor.admr
                rotor_array[i, 27] = curr_rotor.pos.theta_rel(90, "°")
                rotor_array[i, 28] = curr_rotor.pos.theta_rel(180, "°")
                rotor_array[i, 29] = curr_rotor.pos.theta_rel(270, "°")

        return rotor_array
