from .support import Position, Speed
from abc import ABC, abstractmethod
import numpy as np


class BaseDiffuser:

    def __init__(self, main_turbine):

        self.main_turbine = main_turbine
        self.options = self.main_turbine.options
        self.geometry = self.main_turbine.geometry

        self.input_point = self.main_turbine.points[3]
        self.output_point = self.main_turbine.points[4]

        self.static_input_point = self.main_turbine.static_points[3]
        self.static_output_point = self.main_turbine.static_points[4]

        # Tobe updated in self.solve()
        self.m_dot_tot = 0
        self.p_out = 0.

    def solve(self):

        m_ch = self.main_turbine.rotor.m_dot_ch
        n_ch = self.main_turbine.geometry.rotor.n_channels
        __tmp_point = self.input_point.duplicate()

        rho_u_curr = m_ch * n_ch / self.geometry.rotor.a_int
        rho_in = self.input_point.get_variable("rho")
        p_in = self.input_point.get_variable("P")

        dp_max = rho_u_curr ** 2 / rho_in
        p_out_lim = [p_in-dp_max, p_in]

        for i in range(self.options.general_bisection_limit):

            p_out_curr = (p_out_lim[0] + p_out_lim[1]) / 2
            __tmp_point.set_variable("P", p_out_curr)
            __tmp_point.set_variable("s", self.input_point.get_variable("s"))

            h_curr = __tmp_point.get_variable("h")
            rho_curr = __tmp_point.get_variable("rho")

            dh = self.input_point.get_variable("h") - h_curr
            u_enth = np.sqrt(np.abs(2 * dh))
            u_rho = rho_u_curr / rho_curr

            if u_rho > u_enth:

                p_out_lim[1] = p_out_curr

            else:

                p_out_lim[0] = p_out_curr

        p_out_curr = (p_out_lim[0] + p_out_lim[1]) / 2

        self.main_turbine.points[4].set_variable("P", p_out_curr)
        self.main_turbine.points[4].set_variable("s", self.input_point.get_variable("s"))

        self.main_turbine.static_output_point[4].set_variable("P", p_out_curr)
        self.main_turbine.static_output_point[4].set_variable("s", self.input_point.get_variable("s"))


class NoDiffuser(BaseDiffuser):

    """
        Using this class means that the dynamic pressure is considered
        to be completely lost at the outlet of the turbine:

            Total pressure at the outlet of the rotor == static pressure at the outlet of the rotor
            Total enthalpy at the outlet of the rotor is constant

            No residual speed at the outlet: static point = total point

     """

    def solve(self):

        self.main_turbine.points[4].set_variable("P", self.static_input_point.get_variable("P"))
        self.main_turbine.points[4].set_variable("H", self.input_point.get_variable("H"))

        self.main_turbine.points[4].copy_state_to(self.main_turbine.static_points[4])
