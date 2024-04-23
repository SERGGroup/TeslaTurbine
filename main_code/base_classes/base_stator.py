from abc import ABC, abstractmethod
from .support import Position, Speed
from .support import BaseTeslaGeometry
import numpy as np


class BaseStator0D(ABC):

    def __init__(self, main_turbine):

        self.main_turbine = main_turbine
        self.options = self.main_turbine.options.stator
        self.geometry = self.main_turbine.geometry.stator

        self.input_point = self.main_turbine.points[0]
        self.output_point = self.main_turbine.points[1]

        self.static_input_point = self.main_turbine.static_points[0]
        self.static_output_point = self.main_turbine.static_points[1]

        out_pos = Position(self.geometry.r_int, 0)
        self.speed_out = Speed(out_pos)
        self.liq_speed_out = Speed(out_pos)
        self.gas_speed_out = Speed(out_pos)

        self.Ma_1 = 0.
        self.eta_stat = 0.
        self.init_speed = 0.
        self.v_out = 0.
        # self.n = 100
        # self.ax = np.zeros(self.n)
        self.m_dot_s = 0

        self.__m_dot_max = None
        self.__m_dot_max_eval_point = None

        # TODO change in solve
        self.p_out = 0.

    # @property
    # def m_dot_s(self):
    #     return self.m_dot_s

    # @m_dot_s.setter
    # def m_dot_s(self, value):
    #     self.m_dot_s = value

    @abstractmethod
    def solve(self):

        pass

    @property
    def m_dot_max(self):

        if (self.__m_dot_max is None) or not (self.__m_dot_max_eval_point == self.input_point):

            self.__m_dot_max = self.__evaluate_m_dot_max()
            self.__m_dot_max_eval_point = self.input_point.duplicate()

        return self.__m_dot_max

    def __evaluate_m_dot_max(self):

        h_in = self.input_point.get_variable("h")
        s_in = self.input_point.get_variable("s")

        __tmp_point = self.input_point.duplicate()

    def evaluate_two_phase_visc(self, p, x):

        __tmp_point = list()

        for i in range(2):
            __tmp_point.append(self.input_point.duplicate())

        __tmp_point[0].set_variable("p", p)
        __tmp_point[0].set_variable("x", 0)

        __tmp_point[1].set_variable("p", p)
        __tmp_point[1].set_variable("x", 1)

        mu_l = __tmp_point[0].get_variable("visc")
        mu_g = __tmp_point[1].get_variable("visc")

        mu = 0.5 * (mu_l * ((2 * mu_l + mu_g - 2 * x * (mu_l - mu_g)) / (2 * mu_l + mu_g + x * (mu_l - mu_g))) +
                    mu_g * ((2 * mu_g + mu_l - 2 * (1 - x) * (mu_g - mu_l)) / (
                            2 * mu_g + mu_l + (1 - x) * (mu_g - mu_l))))

        return mu, mu_l, mu_g

    def evaluate_lockhart_martinelli(self, p, x):

        epsilon, rho_l, rho_g = self.evaluate_epsilon(p, x)
        mu, mu_l, mu_g = self.evaluate_two_phase_visc(p, x)

        chi_s = ((1 - x) / x) ** 1.8 * (rho_g / rho_l) * (mu_l / mu_g) ** 0.2
        chi = np.sqrt(chi_s)

        return chi, chi_s, mu_l, mu_g, mu

    def evaluate_two_phase_ss(self, p, x):

        __tmp_point = list()

        for i in range(2):
            __tmp_point.append(self.input_point.duplicate())

        __tmp_point[0].set_variable("p", p)
        __tmp_point[0].set_variable("x", 0)

        __tmp_point[1].set_variable("p", p)
        __tmp_point[1].set_variable("x", 1)

        ss_l = __tmp_point[0].get_variable("c")
        ss_g = __tmp_point[1].get_variable("c")

        epsilon, rho_l, rho_g = self.evaluate_epsilon(p, x)

        ss = (epsilon / (ss_g ** 2) + ((1 - epsilon) / ss_l) ** 2 + epsilon *
              (1 - epsilon) * rho_l /(1.35 * p)) ** (-0.5)

        return ss

    def evaluate_epsilon(self, p, x):

        __tmp_point = list()

        for i in range(2):
            __tmp_point.append(self.input_point.duplicate())

        __tmp_point[0].set_variable("p", p)
        __tmp_point[0].set_variable("x", 0)

        __tmp_point[1].set_variable("p", p)
        __tmp_point[1].set_variable("x", 1)

        rho_l = __tmp_point[0].get_variable("rho")
        rho_g = __tmp_point[1].get_variable("rho")

        epsilon = 1 / (1 + (1 - x) / x * (rho_g / rho_l) ** 0.66)

        return epsilon, rho_l, rho_g


class BaseStatorStep(ABC):

    def __init__(self, main_stator, speed: Speed):

        self.main_stator = main_stator
        self.options = main_stator.options

        self.speed = speed
        self.thermo_point = self.main_stator.input_point.duplicate()

    @property
    def pos(self):

        return self.speed.pos

    @property
    def m_dot(self):

        return self.main_stator.m_dot_s

    def get_new_step(self, ds):

        # Evaluate main parameter variation
        dvt, dvr, dp, dh = self.get_variations(ds)

        # Update Position
        new_pos = self.get_new_position(ds)

        # Update Speed
        new_speed = Speed(new_pos)
        vt_new = self.speed.vt + dvt
        vr_new = self.speed.vt + dvr
        new_speed.init_from_codes("vt", vt_new, "vr", vr_new)

        # Init a new step
        new_class = self.self_class()
        new_step = new_class(main_stator=self.main_stator, speed=new_speed)

        # Update Thermodynamic Point
        p_new = self.thermo_point.get_variable("P") + dp
        h_new = self.thermo_point.get_variable("H") + dh
        new_step.thermo_point.set_variable("P", p_new)
        new_step.thermo_point.set_variable("H", h_new)

        return new_step

    @abstractmethod
    def get_new_position(self, ds):

        return self.pos

    @abstractmethod
    def get_variations(self, ds):

        return 0., 0., 0., 0.

    def __init_subclass__(cls, *args, **kwargs):

        def return_self_class(self):

            return type(self)

        cls.self_class = return_self_class


class BaseStator1D(BaseStator0D):

    def __init__(self, main_turbine, stator_step: type(BaseStatorStep)):

        super().__init__(main_turbine)
        self.__stator_step_cls = stator_step
        self.stator_points = list()

        self.__omega = 0.
        self.inlet_point = 0.

    def solve(self):

        # TODO INITIALIZATION
        # flow rate -> initial speed

        first_pos = Position(self.geometry.d_0, omega=0.)
        first_speed = Speed(position=first_pos)
        # Here I have to modify the initialization of first speed, once the mass flow rate is calculated
        # from the previous step
        first_speed.equal_absolute_speed_to(self.main_turbine.stator.speed_out)
        first_step = self.__stator_step_cls(self, first_speed)

        new_step = first_step

        self.stator_points = list()
        ds = (self.geometry.r_0 - self.geometry.r_int) / self.options.n_stator

        for i in range(self.options.n_stator):

            if self.options.profile_stator:
                self.stator_points.append(new_step)

            new_step = new_step.get_new_step(ds)

        if self.options.profile_stator:
            self.stator_points.append(new_step)
