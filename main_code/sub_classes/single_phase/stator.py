from main_code.tesla_turbine_class import TeslaTurbine
from REFPROPConnector import ThermodynamicPoint as TP
import numpy as np

class Stator:

    def __init__(self, main_turbine: TeslaTurbine):

        self.main_turbine = main_turbine
        self.options = self.main_turbine.options.stator

        self.input_point = self.main_turbine.points[0]
        self.output_point = self.main_turbine.points[1]

        self.static_input_point = self.main_turbine.static_points[0]
        self.static_output_point = self.main_turbine.static_points[1]

        self.__tmp_points = list()

        for i in range(3):
            self.__tmp_points.append(self.main_turbine.points[0].duplicate())

        self.geometry = self.main_turbine.geometry.stator

    def __stator_calc(self, n_it):

        phi_n_up = 0.97
        phi_n_down = 0.9

        for i in range(n_it):

            phi_n = (phi_n_up + phi_n_down) / 2

            Ma_1, Xi_diff, Xi_Rodg, Xi_dix, v1, SS_1 = self.__evaluate_phi(phi_n)

            if Xi_diff < 0.0001:
                break
            else:
                if Xi_dix < Xi_Rodg:
                    phi_n_up = phi_n
                else:
                    phi_n_down = phi_n

        return Ma_1, v1, SS_1

    def __evaluate_phi(self, phi_n):

            # [1s] ISOENTROPIC STATOR OUTLET
            self.__tmp_points[0].set_variable("p", self.static_output_point.get_variable("p"))
            self.__tmp_points[0].set_variable("s", self.input_point.get_variable("s"))

            v_1s = np.sqrt(2 * (self.input_point.get_variable("h") - self.__tmp_points[0].get_variable("h")))
            v1 = phi_n * v_1s

            self.output_point.set_variable("h", self.input_point.get_variable("h"))
            self.static_output_point.set_variable("h", self.output_point.get_variable("h") - 0.5 * v1 ** 2)
            T1 = self.static_output_point.get_variable("T")
            x1 = self.static_output_point.get_variable("x")
            rho1 = self.static_output_point.get_variable("rho")

            self.__tmp_points[1].set_variable("p", self.static_output_point.get_variable("p"))
            self.__tmp_points[1].set_variable("x", 0)

            self.__tmp_points[2].set_variable("p", self.static_output_point.get_variable("p"))
            self.__tmp_points[2].set_variable("x", 1)



            mu1 = self.__tmp_points[1].get_variable("visc")
            SS_1 = self.__tmp_points[1].get_variable("c")
            Ma_1 = v1 / SS_1

            # LOSS FACTOR RODGERS
            pitch = 2 * np.pi * self.geometry.r_int / self.geometry.Z_stat
            ni1 = mu1 / rho1
            Re_rodg = v1 * self.geometry.H_s / ni1
            Xi_rodg = (0.05 / (Re_rodg ** 0.2)) * (
                        (3 * np.tan(self.geometry.alpha1) * self.geometry.chord) / pitch + pitch * np.cos(
                    self.geometry.alpha1) / self.geometry.H_s)

            # LOSS FACTOR DIXON (P. 257)

            Xi_dix = 1 / (phi_n ** 2) - 1

            Xi_diff = np.abs(Xi_dix - Xi_rodg)

            return Ma_1, Xi_diff, Xi_dix, Xi_rodg, v1, SS_1


    def solve(self):

        A0 = (2 * np.pi * (1.5 * self.geometry.r_int) / self.geometry.N_s) * self.geometry.Z_stat * self.geometry.H_s
        A1 = self.geometry.Z_stat * self.geometry.throat_width * self.geometry.H_s

        # Point[0] = input_point ----- Total Conditions Inlet Section Data
        self.input_point.set_variable("P", 1)
        self.input_point.set_variable("x", 0)
        h_00 = self.input_point.get_variable("h")
        s_00 = self.input_point.get_variable("s")
        T_00 = self.input_point.get_variable("T")
        rho_00 = self.input_point.get_variable("rho")

        # Static_Point[1] = static_output_point ----- Static Conditions Outlet Section Data
        self.static_output_point.set_variable("P", 1)

        if self.options.iterate_phi:

            Ma_1, v1, SS_1 = self.__stator_calc(n_it=self.options.n_phi_iteration)

        else:

            Ma_1, Xi_diff, Xi_dix, Xi_rodg, v1, SS_1 = self.__evaluate_phi(0.9)

        # Checking Mach Number

        if Ma_1 < 1:
            v_out = v1
        else:
            v_out = SS_1

        # Calculating Static Point [1] Conditions
        self.output_point.set_variable("h", self.input_point.get_variable("h"))
        self.static_output_point.set_variable("h", self.output_point.get_variable("h") - 0.5 * v_out ** 2)

        T1 = self.static_output_point.get_variable("T")
        x1 = self.static_output_point.get_variable("x")
        rho1 = self.static_output_point.get_variable("rho")
        s1 = self.static_output_point.get_variable("s")
        h1 = self.static_output_point.get_variable("h")

        m_dot_s = A1 * rho1 * v_out

         # Calculate Static Point [0] Conditions

        self.static_input_point.set_variable("rho", rho_00)
        v0 = m_dot_s / (rho_00 * A0)
        self.static_input_point.set_variable("h", self.input_point.get_variable("h") - 0.5 * v0 ** 2)

        s0 = self.static_input_point.get_variable("s")
        T0 = self.static_input_point.get_variable("T")
        p0 = self.static_input_point.get_variable("p")
        x0 = self.static_input_point.get_variable("x")
        rho0 = self.static_input_point.get_variable("rho")

          # Calculate Total Point [01] Conditions
        self.output_point.set_variable("s", self.static_output_point.get_variable("s"))
        p_01 = self.output_point.get_variable("p")
        T_01 = self.output_point.get_variable("T")
        rho_01 = self.output_point.get_variable("rho")
        h_01 = self.output_point.get_variable("h")

        Eta_stat = (self.output_point.get_variable("h") - self.static_output_point.get_variable("h"))/(self.output_point.get_variable("h") - self.__tmp_points[0].get_variable("h"))

        return T0, p0, T1, Eta_stat, Ma_1, x0, x1, v1, m_dot_s, h_01, p_01, h1

