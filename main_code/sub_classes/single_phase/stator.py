
from main_code.base_classes import BaseStator0D, BaseStatorStep
import numpy as np


class SPStator(BaseStator0D):


    def __init__(self, main_turbine):

        super().__init__(main_turbine)
        self.__tmp_points = list()

        for i in range(3):
            self.__tmp_points.append(self.main_turbine.points[0].duplicate())

        self.isentropic_output = self.main_turbine.points[0].duplicate()

    def __stator_calc(self, n_it):

        phi_n_up = 0.97
        phi_n_down = 0.9

        for i in range(n_it):

            self.phi_n = (phi_n_up + phi_n_down) / 2

            self.Ma_1, self.Xi_diff, self.Xi_Rodg, self.Xi_dix, v1, SS_1 = self.__evaluate_phi(self.phi_n)

            if self.Xi_diff < 0.0001:
                break
            elif self.Xi_dix > self.Xi_Rodg:
                phi_n_up = self.phi_n
            else:
                phi_n_down = self.phi_n

        return self.Ma_1, v1, SS_1, self.Xi_diff, self.Xi_Rodg, self.Xi_dix

    def __evaluate_phi(self, phi_n):

        self.__tmp_points[0].set_variable("p", self.p_out)
        self.__tmp_points[0].set_variable("s", self.input_point.get_variable("s"))

        dh = self.input_point.get_variable("h") - self.__tmp_points[0].get_variable("h")
        # [1s] ISOENTROPIC STATOR OUTLET
        v_1s = np.sqrt(2 * dh)
        v1 = phi_n * v_1s

        self.output_point.set_variable("h", self.input_point.get_variable("h"))
        self.output_point.set_variable("P", self.p_out)

        self.static_output_point.set_variable("h", self.output_point.get_variable("h") - 0.5 * v1 ** 2)
        self.static_output_point.set_variable("P", self.p_out)

        rho1 = self.static_output_point.get_variable("rho")
        mu1 = self.static_output_point.get_variable("visc")
        SS_1 = self.static_output_point.get_variable("c")
        Ma_1 = v1 / SS_1

        # LOSS FACTOR RODGERS
        pitch = 2 * np.pi * self.geometry.r_int / self.geometry.Z_stat
        ni1 = mu1 / rho1
        Re_rodg = v1 * self.main_turbine.geometry.H_s / ni1
        Xi_rodg = (0.05 / (Re_rodg ** 0.2)) * ((3 * np.tan(self.geometry.alpha_rad) * self.geometry.chord) / pitch + pitch * np.cos(self.geometry.alpha_rad) / self.main_turbine.geometry.H_s)

        # LOSS FACTOR DIXON (P. 257)
        Xi_dix = 1 / (phi_n ** 2) - 1
        Xi_diff = np.abs(Xi_dix - Xi_rodg)

        return Ma_1, Xi_diff, Xi_dix, Xi_rodg, v1, SS_1

    def solve(self):

        self.p_out = self.static_output_point.get_variable("P")

        A0 = (2 * np.pi * (1.5 * self.geometry.r_int) / self.geometry.N_s) * self.geometry.Z_stat * self.main_turbine.geometry.H_s
        A1 = self.geometry.Z_stat * self.main_turbine.geometry.throat_width * self.main_turbine.geometry.H_s

        if self.options.iterate_phi:

            self.Ma_1, v1, SS_1, self.Xi_diff, self.Xi_Rodg, self.Xi_dix, = self.__stator_calc(n_it=self.options.n_phi_iteration)

        else:

            self.Ma_1, self.Xi_diff, self.Xi_dix, self.Xi_rodg, v1, SS_1 = self.__evaluate_phi(0.9)

        # Checking Mach Number

        if self.Ma_1 < 1:
            v_out = v1

        else:
            v_out = SS_1
            self.Ma_1 = 1

        # Calculating Static Point [0] Conditions
        rho_00 = self.input_point.get_variable("rho")
        v0 = self.m_dot_s / (rho_00 * A0)

        self.static_input_point.set_variable("rho", rho_00)
        self.static_input_point.set_variable("h", self.input_point.get_variable("h") - 0.5 * v0 ** 2)

        # Calculating Static Point [1] Conditions
        self.output_point.set_variable("h", self.input_point.get_variable("h"))
        self.output_point.set_variable("P", self.p_out)
        self.static_output_point.set_variable("h", self.output_point.get_variable("h") - 0.5 * v_out ** 2)
        self.static_output_point.set_variable("P", self.p_out)

        rho1 = self.static_output_point.get_variable("rho")
        self.m_dot_s = A1 * rho1 * v_out
        self.speed_out.init_from_codes("v", v_out, "alpha", self.geometry.alpha_rad)

        self.eta_stat = (self.output_point.get_variable("h") - self.static_output_point.get_variable("h"))/(self.output_point.get_variable("h") - self.__tmp_points[0].get_variable("h"))
        self.__tmp_points[0].copy_state_to(self.isentropic_output)

class SPStatorStep(BaseStatorStep):
    def get_new_position(self, ds):
        pass

    def get_variations(self, ds):
        pass

