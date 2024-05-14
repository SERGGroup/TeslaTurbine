from .support import Position, Speed
from abc import ABC, abstractmethod
import numpy as np


class BaseRotorStep(ABC):

    def __init__(self, main_rotor, speed: Speed):

        self.main_rotor = main_rotor
        self.options = main_rotor.options
        self.geometry = self.main_rotor.geometry

        self.speed = speed
        self.total_point = self.main_rotor.input_point.duplicate()
        self.static_point = self.main_rotor.input_point.duplicate()

        self.rothalphy = self.main_rotor.inlet_rothalpy
        self.__omega = 0.

    @property
    def pos(self):

        return self.speed.pos

    @property
    def m_dot(self):

        return self.main_rotor.m_dot_ch

    def get_new_step(self, dr):

        # Update Position
        self.new_pos = self.speed.get_new_position(dr)
        dvr, dvt, dp, dh = self.get_variations(dr)

        # Update Speed
        new_speed = Speed(self.new_pos)
        vt_new = self.speed.vt + dvt
        vr_new = self.speed.vr + dvr
        new_speed.init_from_codes("vt", vt_new, "vr", vr_new)

        # Init a new step
        new_class = self.self_class()
        new_step = new_class(main_rotor=self.main_rotor, speed=new_speed)

        # Update Enthalpy through Rothalpy Conservation
        new_step.rothalpy = self.rothalphy + dh
        h_new = new_step.rothalpy + new_speed.u * new_speed.vt - 1/2 * new_speed.v ** 2

        # Update Thermodynamic Point
        p_new = self.static_point.get_variable("P") + dp
        new_step.static_point.set_variable("P", p_new)
        new_step.static_point.set_variable("H", h_new)

        if self.options.evaluate_internal_total_condition:
            total_point = new_step.static_point.get_stagnation_point(speed=new_speed.v)
            total_point.copy_state_to(new_step.total_point)

        return new_step

    @abstractmethod
    def get_variations(self, dr):

        """Evaluate changes in speed, STATIC pressure and TOTAL enthalpy (for adiabatic expansion consider dh = 0)"""
        return 0., 0., 0., 0.

    def __init_subclass__(cls, *args, **kwargs):

        def return_self_class(self):

            return type(self)

        cls.self_class = return_self_class


class BaseRotor(ABC):

    __dv_perc = 0.1315
    __omega = None

    def __init__(

            self,
            main_turbine,
            rotor_step: type(BaseRotorStep)

    ):

        self.eta_tt = 0.
        self.eta_ss = 0.
        self.inlet_rothalpy = 0.
        self.gap_losses_control = True

        self.main_turbine = main_turbine
        self.geometry = self.main_turbine.geometry.rotor
        self.options = self.main_turbine.options.rotor

        self.input_point = self.main_turbine.points[2]
        self.output_point = self.main_turbine.points[3]

        self.rotor_points = list()
        self.rotor_step_cls = rotor_step

        inlet_pos = Position(self.geometry.r_out, 0)
        self.rotor_inlet_speed = Speed(inlet_pos)
        self.first_speed = Speed(inlet_pos)

    def solve(self):

        self.rotor_points = list()
        self.evaluate_gap_losses()

        first_pos = Position(self.geometry.r_out, self.omega)
        self.first_speed = Speed(position=first_pos)
        self.first_speed.equal_absolute_speed_to(self.rotor_inlet_speed)
        self.inlet_rothalpy = self.main_turbine.points[2].get_variable("h") - self.first_speed.vt * self.first_speed.u

        first_step = self.rotor_step_cls(self, self.first_speed)
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

        self.rotor_points[-1].static_point.copy_state_to(self.main_turbine.static_points[3])

        if self.options.evaluate_internal_total_condition:
            self.rotor_points[-1].total_point.copy_state_to(self.main_turbine.points[3])

        else:
            new_speed = self.rotor_points[-1].speed.v
            total_point = self.rotor_points[-1].static_point.get_stagnation_point(speed=new_speed)
            total_point.copy_state_to(self.main_turbine.points[3])

    @abstractmethod
    def evaluate_gap_losses(self):

        """

            This Function must update:
                - self.rotor_inlet_speed
                - self.main_turbine.points[2]

        """

        # No Loss Model (The total properties will be preserved)
        self.main_turbine.points[1].copy_state_to(self.main_turbine.points[2])

        self.rotor_inlet_speed = self.main_turbine.stator.speed_out
        rho_out = self.main_turbine.points[1].get_variable("rho")

        self.main_turbine.points[2].set_variable("P", 101325)
        self.main_turbine.points[2].set_variable("h", 500000)

    def evaluate_performances(self):

        self.isentropic_outlet = self.static_points[2].duplicate()

        self.isentropic_outlet.set_variable("s", self.static_points[2].get_variable("s"))
        self.isentropic_outlet.set_variable("p", self.static_points[3].get_variable("p"))

        self.Eta_rotor_ss = (self.static_points[2].get_variable("h") - self.static_points[3].get_variable("h")) / (
                self.static_points[2].get_variable("h") - self.isentropic_outlet.get_variable("h"))

        self.work = self.rotor.first_speed.vt * self.rotor.first_speed.u - self.rotor.rotor_points[-1].speed.vt * \
                    self.rotor.rotor_points[-1].speed.u
        self.power = self.work * self.rotor.m_dot_ch

    @property
    def m_dot_ch(self):

        return self.main_turbine.stator.m_dot_s / self.geometry.n_channels

    def get_rotor_array(self):

        rotor_array = np.empty((len(self.rotor_points), 30))
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

                rotor_array[i, 11] = curr_rotor.total_point.get_variable("P")
                rotor_array[i, 12] = curr_rotor.total_point.get_variable("h")
                rotor_array[i, 13] = curr_rotor.total_point.get_variable("T")
                rotor_array[i, 14] = curr_rotor.total_point.get_variable("s")
                rotor_array[i, 15] = curr_rotor.total_point.get_variable("rho")
                rotor_array[i, 16] = curr_rotor.total_point.get_variable("mu")

                # rotor_array[i, 17] = curr_rotor.Ma_R_a
                # rotor_array[i, 18] = curr_rotor.Ma_R_r

                # rotor_array[i, 19] = curr_rotor.Re_a
                rotor_array[i, 20] = curr_rotor.total_point.get_variable("quality")
                # rotor_array[i, 21] = curr_rotor.cosy

                rotor_array[i, 24] = curr_rotor.pos.theta_rel(0, "°")
                rotor_array[i, 25] = curr_rotor.pos.gamma_rel(0, "°")
                # rotor_array[i, 26] = curr_rotor.admr
                rotor_array[i, 27] = curr_rotor.pos.theta_rel(90, "°")
                rotor_array[i, 28] = curr_rotor.pos.theta_rel(180, "°")
                rotor_array[i, 29] = curr_rotor.pos.theta_rel(270, "°")

        return rotor_array

    @property
    def dv_perc(self):

        if self.__omega is None:

            return self.__dv_perc

        else:

            try:

                return self.rotor_inlet_speed.vt / self.omega - self.geometry.r_out - 1

            except:

                return None

    @dv_perc.setter
    def dv_perc(self, dv_perc):

        self.__dv_perc = dv_perc
        self.__omega = None

    @property
    def omega(self):

        if self.__omega is None:

            try:

                return self.rotor_inlet_speed.vt / ((self.dv_perc + 1) * self.geometry.r_out)

            except:

                return None

        else:

            return self.__omega

    @omega.setter
    def omega(self, omega):

        self.__omega = omega

    @property
    def rpm(self):

        return self.omega * 60 / (2 * np.pi)

    @rpm.setter
    def rpm(self, rpm):

        self.omega = rpm / 60 * (2 * np.pi)