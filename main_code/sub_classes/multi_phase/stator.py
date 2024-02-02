#from main_code.tesla_turbine_class import TeslaTurbine
from REFPROPConnector import ThermodynamicPoint as TP
from main_code.base_classes.support import Speed, Position
import numpy as np


class Stator:

    def __init__(self, main_turbine):

        self.main_turbine = main_turbine
        self.options = self.main_turbine.options.stator

        self.input_point = self.main_turbine.points[0]
        self.output_point = self.main_turbine.points[1]

        self.static_input_point = self.main_turbine.static_points[0]
        self.static_output_point = self.main_turbine.static_points[1]
        self.p_out = self.static_output_point.get_variable("P")

        self.__tmp_points = list()

        for i in range(3):
            self.__tmp_points.append(self.main_turbine.points[0].duplicate())

        self.geometry = self.main_turbine.geometry.stator
        pos = Position(self.geometry.r_int, 0)
        self.speed_out = Speed(pos)

        self.Ma_1 = 0.
        self.m_dot_s = 0.
        self.eta_stat = 0.
        self.phi_n = 0.9

    def __stator_calc(self, n_it):

        phi_n_up = 0.97
        phi_n_down = 0.9

        for i in range(n_it):

            self.phi_n = (phi_n_up + phi_n_down) / 2

            Ma_1, Xi_diff, Xi_Rodg, Xi_dix, v1, SS_1 = self.__evaluate_phi(self.phi_n)

            if Xi_diff < 0.0001:
                break
            else:
                 if Xi_dix < Xi_Rodg:
                    phi_n_down = self.phi_n
                 else:
                    phi_n_up = self.phi_n

        return Ma_1, v1, SS_1

    def __evaluate_phi(self, phi_n):

        self.__tmp_points[0].set_variable("p", self.static_output_point.get_variable("p"))
        self.__tmp_points[0].set_variable("s", self.input_point.get_variable("s"))

        # [1s] ISOENTROPIC STATOR OUTLET
        dh = self.input_point.get_variable("h") - self.__tmp_points[0].get_variable("h")
        v_1s = np.sqrt(2 * dh)
        v1 = phi_n * v_1s

        self.output_point.set_variable("h", self.input_point.get_variable("h"))
        self.output_point.set_variable("P", self.p_out)

        self.static_output_point.set_variable("h", self.output_point.get_variable("h") - 0.5 * v1 ** 2 )
        self.static_output_point.set_variable("P", self.p_out)

        x1 = self.static_output_point.get_variable("x")
        rho1 = self.static_output_point.get_variable("rho")

        self.__tmp_points[1].set_variable("p", self.static_output_point.get_variable("p"))
        self.__tmp_points[1].set_variable("x", 0)

        self.__tmp_points[2].set_variable("p", self.static_output_point.get_variable("p"))
        self.__tmp_points[2].set_variable("x", 1)

        mu_1l = self.__tmp_points[1].get_variable("visc")
        mu_1g = self.__tmp_points[2].get_variable("visc")

        mu1 = 0.5 * (mu_1l * ((2 * mu_1l + mu_1g - 2 * x1 * (mu_1l - mu_1g))/(2 * mu_1l + mu_1g + x1 * (mu_1l - mu_1g)) + mu_1g * ((2 * mu_1g + mu_1l - 2 * (1 - x1)*(mu_1g - mu_1l))/(2 * mu_1g + mu_1l + (1 - x1) * (mu_1g - mu_1l)))))

        rho_1l = self.__tmp_points[1].get_variable("rho")
        rho_1g = self.__tmp_points[2].get_variable("rho")

        void_f1 = 1/(1 + ((1-x1)/x1) * (rho_1g/rho_1l) ** 0.66)

        SS_1l = self.__tmp_points[1].get_variable("c")
        SS_1g = self.__tmp_points[2].get_variable("c")

        SS_1 = (void_f1/(SS_1g ** 2) + (1-void_f1 ** 2)/(SS_1l ** 2) + void_f1 * (1 - void_f1) * rho_1l/(1.35 * self.static_output_point.get_variable("p"))) ** (-0.5)
        Ma_1 = v1/SS_1

        # LOSS FACTOR RODGERS
        pitch = 2 * np.pi * self.geometry.r_int / self.geometry.Z_stat
        ni1 = mu1 / rho1
        Re_rodg = v1 * self.geometry.H_s / ni1
        Xi_rodg = (0.05/(Re_rodg ** 0.2)) * ((3 * np.tan(self.geometry.alpha_rad) * self.geometry.chord) / pitch + pitch * np.cos(self.geometry.alpha_rad) / self.geometry.H_s)

        # LOSS FACTOR DIXON (P. 257)

        Xi_dix = 1 / (phi_n ** 2) - 1

        Xi_diff = np.abs(Xi_dix - Xi_rodg)

        return Ma_1, Xi_diff, Xi_dix, Xi_rodg, v1, SS_1

    def solve(self):

        self.p_out = self.static_output_point.get_variable("P")

        A0 = (2 * np.pi * (1.5 * self.geometry.r_int) / self.geometry.N_s) * self.geometry.Z_stat * self.geometry.H_s
        A1 = self.geometry.Z_stat * self.geometry.throat_width * self.geometry.H_s

        if self.options.iterate_phi:

            self.Ma_1, v1, SS_1 = self.__stator_calc(n_it=self.options.n_phi_iteration)

        else:

            self.Ma_1, Xi_diff, Xi_dix, Xi_rodg, v1, SS_1 = self.__evaluate_phi(self.phi_n)

        # Checking Mach Number

        if self.Ma_1 < 1:
            v_out = v1
        else:
            v_out = SS_1
            self.Ma_1 = 1


        # Calculating Static Point [1] Conditions
        self.output_point.set_variable("h", self.input_point.get_variable("h"))
        self.output_point.set_variable("P", self.p_out)

        self.static_output_point.set_variable("h", self.output_point.get_variable("h") - 0.5 * v_out ** 2)
        self.static_output_point.set_variable("P", self.p_out)


        x1 = self.static_output_point.get_variable("x")
        rho1 = self.static_output_point.get_variable("rho")

        mu_1l = self.__tmp_points[1].get_variable("visc")
        mu_1g = self.__tmp_points[2].get_variable("visc")

        mu1 = 0.5 * (mu_1l * ((2 * mu_1l + mu_1g - 2 * x1 * (mu_1l - mu_1g))/(2 * mu_1l + mu_1g + x1 * (mu_1l - mu_1g)) + mu_1g * ((2 * mu_1g + mu_1l - 2 * (1 - x1)*(mu_1g - mu_1l))/(2 * mu_1g + mu_1l + (1 - x1) * (mu_1g - mu_1l)))))

        rho_1l = self.__tmp_points[1].get_variable("rho")
        rho_1g = self.__tmp_points[2].get_variable("rho")

        void_f1 = 1/(1 + ((1-x1)/x1) * (rho_1g/rho_1l) ** 0.66)

        SS_1l = self.__tmp_points[1].get_variable("c")
        SS_1g = self.__tmp_points[2].get_variable("c")

        SS_1 = (void_f1/(SS_1g ** 2) + (1-void_f1 ** 2)/(SS_1l ** 2) + void_f1 * (1 - void_f1) * rho_1l/(1.35 * self.static_output_point.get_variable("p"))) ** (-0.5)
        Ma_1 = v_out/SS_1

        self.m_dot_s = A1 * rho1 * v_out

        # Calculating Static Point [0] Conditions
        rho_00 = self.input_point.get_variable("rho")
        v0 = self.m_dot_s / (rho_00 * A0)

        self.static_input_point.set_variable("rho", rho_00)
        self.static_input_point.set_variable("h", self.input_point.get_variable("h") - 0.5 * v0 ** 2)

        self.speed_out = Speed(Position(self.geometry.r_int, 0))
        self.speed_out.init_from_codes("v", v_out, "alpha", self.geometry.alpha_rad)

        dh_st = self.output_point.get_variable("h") - self.static_output_point.get_variable("h")
        dh_is = self.output_point.get_variable("h") - self.__tmp_points[0].get_variable("h")
        self.eta_stat = dh_st/dh_is






