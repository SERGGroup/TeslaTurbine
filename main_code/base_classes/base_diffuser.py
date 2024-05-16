from .support import Position, Speed
from abc import ABC, abstractmethod
import numpy as np


class BaseDiffuser:

    """
        Using this class means that the dynamic pressure at the outlet of the rotor is considered
        to be completely lost. Moreover, the outlet speed is evaluated from the internal radius dimension:

            - Total pressure in the intermediate point (the outlet of the rotor) == static pressure
            - Total enthalpy at the outlet of the rotor is constant (the compression energy is converted in vorticity)

            - No residual in the intermediate point: static point = total point

            - from the intermediate point to the outlet the total pressure and enthalpy are kept constant (no losses)
            - The static points depends on the outlet density (has to be iterated)

     """

    def __init__(self, main_turbine):

        self.main_turbine = main_turbine
        self.options = self.main_turbine.options.diffuser
        self.geometry = self.main_turbine.geometry

        self.input_point = self.main_turbine.points[3]
        self.output_point = self.main_turbine.points[4]

        self.static_input_point = self.main_turbine.static_points[3]
        self.static_output_point = self.main_turbine.static_points[4]
        self.intermediate_point = self.main_turbine.points[3].duplicate()

        # Tobe updated in self.solve()
        self.m_dot_tot = 0
        self.p_out = 0.

    def solve(self):

        # Total pressure in the intermediate point (the outlet of the rotor) == static pressure
        # Total enthalpy at the outlet of the rotor is constant (the compression energy is converted in vorticity)
        # No residual in the intermediate point: static point = total point
        self.intermediate_point.set_variable("P", self.static_input_point.get_variable("P"))
        self.intermediate_point.set_variable("H", self.input_point.get_variable("H"))

        # From the intermediate point to the outlet the total pressure and enthalpy are kept constant (no losses)
        self.intermediate_point.copy_state_to(self.main_turbine.points[4])

        # The static points depends on the outlet density (has to be iterated)
        rho_u_curr = self.main_turbine.m_dot_tot / self.geometry.rotor.a_int
        rho_0 = self.intermediate_point.get_variable("rho")
        __tmp_point = self.intermediate_point.duplicate()

        rho_curr = rho_0
        for i in range(self.options.n_max_iter):

            u_curr = rho_u_curr / rho_curr
            __tmp_point = self.intermediate_point.get_static_point(speed=u_curr)

            rho_old = rho_curr
            rho_curr = __tmp_point.get_variable("rho")
            err = np.abs(rho_curr - rho_old) / rho_old

            if err < self.options.tol:

                break

            else:

                rho_curr = rho_old * (1 - self.options.alpha) + self.options.alpha * __tmp_point.get_variable("rho")

        __tmp_point.copy_state_to(self.main_turbine.static_output_point[4])


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
