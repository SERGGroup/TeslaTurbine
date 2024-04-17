
from main_code.base_classes import BaseStator0D, BaseStatorStep
import numpy as np


class SPStator(BaseStator0D):


    def __init__(self, main_turbine):

        super().__init__(main_turbine)
        self.__tmp_points = list()

        for i in range(3):
            self.__tmp_points.append(self.main_turbine.points[0].duplicate())

        self.isentropic_output = self.main_turbine.points[0].duplicate()
        self.v_1s = 0.

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
        self.v_1s = np.sqrt(2 * dh)
        v1 = phi_n * self.v_1s

        self.static_output_point.set_variable("h", self.input_point.get_variable("h") - 0.5 * v1 ** 2)
        self.static_output_point.set_variable("P", self.p_out)

        self.output_point.set_variable("h", self.input_point.get_variable("h"))
        self.output_point.set_variable("P", self.p_out + 0.5 * self.static_output_point.get_variable("rho") * v1 ** 2)


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
            self.v_out = v1

        else:
            self.v_out = SS_1
            self.Ma_1 = 1

        # Calculating Static Point [0] Conditions
        rho_00 = self.input_point.get_variable("rho")
        v0 = self.m_dot_s / (rho_00 * A0)

        self.static_input_point.set_variable("rho", rho_00)
        self.static_input_point.set_variable("h", self.input_point.get_variable("h") - 0.5 * v0 ** 2)

        # Calculating Static Point [1] Conditions

        self.static_output_point.set_variable("h", self.input_point.get_variable("h") - 0.5 * self.v_out ** 2)
        self.static_output_point.set_variable("P", self.p_out)

        self.output_point.set_variable("h", self.input_point.get_variable("h"))
        self.output_point.set_variable("P", self.p_out + 0.5 * self.static_output_point.get_variable("rho") * self.v_out ** 2)

        rho1 = self.static_output_point.get_variable("rho")
        self.m_dot_s = A1 * rho1 * self.v_out
        self.speed_out.init_from_codes("v", self.v_out, "alpha", self.geometry.alpha_rad)

        self.eta_stat = (self.output_point.get_variable("h") - self.static_output_point.get_variable("h"))/(self.output_point.get_variable("h") - self.__tmp_points[0].get_variable("h"))
        self.__tmp_points[0].copy_state_to(self.isentropic_output)


class SPStatorMil(BaseStator0D):

    def __init__(self, main_turbine):

        super().__init__(main_turbine)

        self.stator_eff = 0.
        self.stator_mil_in = self.main_turbine.points[0].duplicate()
        self.stator_mil_out = self.main_turbine.points[0].duplicate()
        self.isentropic_output = self.main_turbine.points[0].duplicate()
        self.n = 100

        self.phi_n = 0.
        self.v1_ss = 0.
        self.out_speed = 0.

    def initialization(self):

        self.p_out = self.main_turbine.static_points[1].get_variable("P")

        self.stator_mil_in.set_variable("T", self.main_turbine.T_in)
        self.stator_mil_in.set_variable("P", self.main_turbine.P_in)

        self.isentropic_output.set_variable("p", self.p_out)
        self.isentropic_output.set_variable("s", self.stator_mil_in.get_variable("s"))

        self.v1_ss = np.sqrt(2 * (self.stator_mil_in.get_variable("h") - self.isentropic_output.get_variable("h")))

        h_out = self.stator_mil_in.get_variable("h") - self.stator_eff * (self.stator_mil_in.get_variable("h")
                                                                          - self.isentropic_output.get_variable("h"))
        self.stator_mil_out.set_variable("P", self.p_out)
        self.stator_mil_out.set_variable("h", h_out)

    def sonic_cond_evaluation(self, n):

        self.initialization()

        Hh = np.linspace(self.stator_mil_in.get_variable("h"), self.stator_mil_out.get_variable("h"), n)
        Ph = np.linspace(self.stator_mil_in.get_variable("p") - 100, self.stator_mil_out.get_variable("p"), n)
        Ss = np.linspace(self.stator_mil_in.get_variable("s"), self.stator_mil_in.get_variable("s"), n)

        A_throat = self.main_turbine.geometry.throat_width * self.main_turbine.geometry.H_s

        rhov_throat = 0.
        i_max = 0.
        v_sound = 0.

        v = np.zeros(n)
        rhov = np.zeros(n)

        m_dot_arr = np.zeros(n)

        int_points = list()

        for i in range(n):
            int_points.append(self.stator_mil_in.duplicate())

        for i in range(n):
            int_points[i].set_variable("P", Ph[i])
            int_points[i].set_variable("s", Ss[i])

            h_curr = (self.stator_mil_in.get_variable("h") - self.stator_eff *
                      (self.stator_mil_in.get_variable("h") - int_points[i].get_variable("h")))

            int_points[i].set_variable("P", Ph[i])
            int_points[i].set_variable("h", h_curr)

            v[i] = np.sqrt(2 * (self.stator_mil_in.get_variable("h") - h_curr))
            rhov[i] = v[i] * int_points[i].get_variable("rho")

            if rhov[i] > rhov_throat:

                i_max = i
                v_sound = v[i]
                rhov_throat = rhov[i]

            m_dot_s_in = self.main_turbine.geometry.stator.Z_stat * A_throat * rhov_throat
            m_dot_arr[i] = m_dot_s_in

        h_out = self.stator_mil_in.get_variable("h") - (v_sound ** 2) / 2
        v = v_sound

        self.m_dot_s = self.main_turbine.geometry.stator.Z_stat * A_throat * rhov_throat

        return v_sound, h_out

    def solve(self):

        self.out_speed, h1 = self.sonic_cond_evaluation(self.n)

        self.phi_n = self.out_speed / self.v1_ss
        self.Ma_1 = self.out_speed / self.stator_mil_out.get_variable("c")

        self.stator_mil_in.copy_state_to(self.main_turbine.points[0])
        self.stator_mil_in.copy_state_to(self.main_turbine.static_points[0])

        self.main_turbine.static_points[1].set_variable("h", h1)
        self.main_turbine.static_points[1].set_variable("p", self.p_out)

        self.speed_out.init_from_codes("v", self.out_speed, "alpha", self.geometry.alpha_rad)

        self.main_turbine.points[1].set_variable("h", self.main_turbine.points[0].get_variable("h"))
        # self.main_turbine.points[1].set_variable("s", self.main_turbine.static_points[1].get_variable("s"))
        self.main_turbine.points[1].set_variable("p", self.main_turbine.static_points[1].get_variable("p") + 0.5 *
                                                 self.main_turbine.static_points[1].get_variable("rho") * self.out_speed
                                                 ** 2)


class SPStatorStep(BaseStatorStep):
    def get_new_position(self, ds):
        pass

    def get_variations(self, ds):
        pass

