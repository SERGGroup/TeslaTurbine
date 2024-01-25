from code.tesla_turbine_class import TeslaTurbine
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

         # p_try = self.output_point.get_variable("P")

        phi_n_up = 0.97
        phi_n_down = 0.9

        for i in range(n_it):

            phi_n = (phi_n_up + phi_n_down) / 2

            Ma_1, Xi_diff = self.__evaluate_phi(phi_n)

            # Update

        return Ma_1

    def __evaluate_phi(self, phi_n):

        # [1s] ISOENTROPIC STATOR OUTLET
        self.__tmp_points[0].set_variable("p", self.static_output_point.get_variable("p"))
        self.__tmp_points[0].set_variable("s", self.input_point.get_variable("s"))

        v_1s = np.sqrt(2 * (self.input_point.get_variable("h") - self.__tmp_points[0].get_variable("h")))
        v1 = phi_n * v_1s

        self.output_point.set_variable("h", self.input_point.get_variable("h"))
        self.static_output_point.set_variable("h", self.output_point.get_variable("h") - 0.5 * v1 ** 2 )
        T1 = self.static_output_point.get_variable("T")
        x1 = self.static_output_point.get_variable("x")

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

        Xi_diff = 1

        return Ma_1, Xi_diff


    def solve(self):

        if self.options.iterate_phi:

            Ma_1 = self.__stator_calc(n_it=self.options.n_phi_iteration)

        else:

            Ma_1, Xi_diff = self.__evaluate_phi(0.9)


