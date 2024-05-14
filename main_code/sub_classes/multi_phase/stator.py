from main_code.base_classes import BaseStator0D, BaseStatorStep, BaseStator1D
from main_code.base_classes.support import Speed, Position
import numpy as np


class TPStator0D(BaseStator0D):

    phi_n = 0.9

    def __init__(self, main_turbine):

        super().__init__(main_turbine)
        self.__tmp_points = list()


        for i in range(3):
            self.__tmp_points.append(self.main_turbine.points[0].duplicate())

        self.v_out = 0.

    def __stator_calc(self, n_it):

        phi_n_up = 0.97
        phi_n_down = 0.9
        Ma_1 = 0.
        v1 = 0.
        SS_1 = 0.

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

        mu1, mu_l, mu_g = self.evaluate_two_phase_visc(self.static_output_point.get_variable("p"),
                                                       self.static_output_point.get_variable("x"))
        ss1 = self.evaluate_two_phase_ss(self.static_output_point.get_variable("p"),
                                         self.static_output_point.get_variable("x"))
        Ma_1 = v1/ss1

        # LOSS FACTOR RODGERS
        pitch = 2 * np.pi * self.geometry.r_int / self.geometry.Z_stat
        ni1 = mu1 / self.static_output_point.get_variable("rho")
        Re_rodg = v1 * self.geometry.H_s / ni1
        Xi_rodg = (0.05/(Re_rodg ** 0.2)) * ((3 * np.tan(self.geometry.alpha_rad) * self.geometry.chord) / pitch + pitch * np.cos(self.geometry.alpha_rad) / self.geometry.H_s)

        # LOSS FACTOR DIXON (P. 257)

        Xi_dix = 1 / (phi_n ** 2) - 1

        Xi_diff = np.abs(Xi_dix - Xi_rodg)

        return Ma_1, Xi_diff, Xi_dix, Xi_rodg, v1, ss1

    def solve(self):

        self.p_out = self.static_output_point.get_variable("P")

        A0 = (2 * np.pi * (1.5 * self.geometry.r_int) / self.geometry.N_s) * self.geometry.Z_stat * self.geometry.H_s
        A1 = self.geometry.Z_stat * self.main_turbine.geometry.throat_width * self.geometry.H_s

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

        mu1 = self.evaluate_two_phase_visc(self.static_output_point.get_variable("p"),
                                             self.static_output_point.get_variable("x"))
        ss1 = self.evaluate_two_phase_ss(self.static_output_point.get_variable("p"),
                                           self.static_output_point.get_variable("x"))
        self.Ma_1 = v1/ss1

        self.m_dot_s = A1 * self.static_output_point.get_variable("rho") * v_out

        # Calculating Static Point [0] Conditions
        rho_00 = self.input_point.get_variable("rho")
        v0 = self.m_dot_s/ (rho_00 * A0)

        self.static_input_point.set_variable("rho", rho_00)
        self.static_input_point.set_variable("h", self.input_point.get_variable("h") - 0.5 * v0 ** 2)

        self.speed_out = Speed(Position(self.geometry.r_int, 0))
        self.speed_out.init_from_codes("v", v_out, "alpha", self.geometry.alpha_rad)

        dh_st = self.output_point.get_variable("h") - self.static_output_point.get_variable("h")
        dh_is = self.output_point.get_variable("h") - self.__tmp_points[0].get_variable("h")
        self.eta_stat = dh_st/dh_is

        self.p_out = self.output_point.get_variable("p")
        self.v_out = v_out
        self.H_1 = self.static_output_point.get_variable("h")


class TPStatorMil(BaseStator0D):

    def __init__(self, main_turbine):

        super().__init__(main_turbine)

        self.stator_eff = 0.
        self.stator_mil_in = self.main_turbine.points[0].duplicate()
        self.stator_mil_out = self.main_turbine.points[0].duplicate()
        self.isentropic_output = self.main_turbine.points[0].duplicate()
        self.n = 1000

        self.phi_n = 0.
        self.v1_ss = 0.
        self.x_in = 0.
        self.x_out = 0.

    def initialization(self):

        self.p_out = self.main_turbine.static_points[1].get_variable("P")

        if self.main_turbine.points[0].get_variable("x") < 0 or self.main_turbine.points[0].get_variable("x") > 1:
            self.stator_mil_in.set_variable("T", self.main_turbine.points[0].get_variable("T"))
        else:
            self.stator_mil_in.set_variable("x", self.main_turbine.points[0].get_variable("x"))

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

        # Hh = np.linspace(self.stator_mil_in.get_variable("h"), self.stator_mil_out.get_variable("h"), n)
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

            if self.main_turbine.options.stator.metastability_check is True:
                int_points[i].metastability = "liq"
            else:
                int_points[i].metastability = "Equilibrium"

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
        self.m_dot_s = self.main_turbine.geometry.stator.Z_stat * A_throat * rhov_throat

        return v_sound, h_out

    def solve(self):

        out_speed, h1 = self.sonic_cond_evaluation(self.n)

        self.phi_n = out_speed / self.v1_ss

        self.main_turbine.static_points[1].set_variable("h", h1)
        self.main_turbine.static_points[1].set_variable("p", self.p_out)
        total_point = self.main_turbine.static_points[1].get_stagnation_point(out_speed)
        total_point.copy_state_to(self.main_turbine.points[1])

        self.speed_out.init_from_codes("v", out_speed, "alpha", self.geometry.alpha_rad)

        if self.main_turbine.points[0].get_variable("x") > 0:
            self.x_in = self.main_turbine.points[0].get_variable("x")
        else:
            self.x_in = 0
        self.x_out = self.main_turbine.points[1].get_variable("x")

        # SPEED CALCULATIONS for Inlet Rotor Requirements
        self.speed_out = Speed(Position(self.geometry.r_int, 0))
        self.speed_out.init_from_codes("v", out_speed, "alpha", self.geometry.alpha_rad)

        # STATIC STATOR INLET CONDITIONS for Performance Evaluation
        v0 = self.m_dot_s / (self.main_turbine.points[0].get_variable("rho") * 2 * np.pi * self.geometry.r_0 * self.main_turbine.geometry.H_s * self.geometry.Z_stat)
        static_point = self.main_turbine.points[0].get_static_point(v0)
        static_point.copy_state_to(self.main_turbine.static_points[0])


class TPStatorStep(BaseStatorStep):
    def __init__(self, main_stator, speed: Speed):

        super().__init__(main_stator, speed)

        self.liq_speed = Speed(speed.pos)
        self.vap_speed = Speed(speed.pos)

        self.liq_phase = self.main_stator.input_point.duplicate()
        self.vap_phase = self.main_stator.input_point.duplicate()

    def get_new_position(self, ds):
        pass

    def get_variations(self, ds):
        pass


class TPStator1D(BaseStator1D):

    def __init__(self, main_turbine):

        super().__init__(main_turbine, BaseStatorStep)

        self.main_turbine = main_turbine
        self.options = self.main_turbine.options.stator

        self.input_point = self.main_turbine.points[0]
        self.output_point = self.main_turbine.points[1]

        self.static_input_point = self.main_turbine.static_points[0]
        self.static_output_point = self.main_turbine.static_points[1]

        self.p_out = self.static_output_point.get_variable("P")

        self.geometry = self.main_turbine.geometry.stator
        pos = Position(self.geometry.r_int, 0)
        self.speed_out = Speed(pos)

        self.IT = 10
        self.cosyy = np.zeros(self.IT)
        self.m_dot_s = 0.

        self.h = np.zeros(self.options.n_stator)
        self.s = np.zeros(self.options.n_stator)
        self.rho = np.zeros(self.options.n_stator)
        self.rho_l = np.zeros(self.options.n_stator)
        self.rho_g = np.zeros(self.options.n_stator)
        self.x = np.zeros(self.options.n_stator)
        self.p = np.zeros(self.options.n_stator)
        self.epsilon = np.zeros(self.options.n_stator)
        self.ss = np.zeros(self.options.n_stator)
        self.Ma = np.zeros(self.options.n_stator)
        self.mu = np.zeros(self.options.n_stator)
        self.mu_l = np.zeros(self.options.n_stator)
        self.mu_g = np.zeros(self.options.n_stator)
        self.chi = np.zeros(self.options.n_stator)
        self.chi_s = np.zeros(self.options.n_stator)

        self.cosy = np.zeros(self.options.n_stator)
        self.phi_g = np.zeros(self.options.n_stator)
        self.phi_l = np.zeros(self.options.n_stator)
        self.c1 = np.zeros(self.options.n_stator)
        self.c2 = np.zeros(self.options.n_stator)
        self.dh_l = np.zeros(self.options.n_stator)
        self.dh_g = np.zeros(self.options.n_stator)
        self.Re_l = np.zeros(self.options.n_stator)
        self.Re_g = np.zeros(self.options.n_stator)

        self.dpdl_l = np.zeros(self.options.n_stator)
        self.dpdl_fr = np.zeros(self.options.n_stator)
        self.f_l = np.zeros(self.options.n_stator)

        self.__stator_points = list()

        for i in range(self.options.n_stator):

            self.__stator_points.append(self.main_turbine.points[0].duplicate())

    def churchill_friction_factor(self, Re, rough, D):

        A = (2.457 * np.log(1 / ((7 / Re) ** 0.9 + 0.27 * rough / D))) ** 16
        B = (37530 / Re) ** 16

        f = 2 * ((8 / Re) ** 12 + 1 / ((A + B) ** 1.5)) ** (1 / 12)

        return f

    def solve(self):

        pass

    def solve_from_flow_rate(self):

        x_ss, y_ss, y_ps, x_ps, m, alpha, alpha_vert, dl = self.geometry.get_coordinates(self.options.n_stator)
        A_th, A_eff = self.geometry.get_area(x_ss, y_ss, y_ps, x_ps, m, self.options.n_stator)

        # Calculating Velocity and Thermodynamic Static Conditions at Stator Inlet [0]

        rho_try = self.input_point.get_variable("rho")
        first_speed = self.m_dot_s / (self.input_point.get_variable("rho") * A_eff[0])
        static_point = self.input_point.get_static_point(first_speed)
        static_point.copy_state_to(self.main_turbine.static_points[0])

        if self.main_turbine.static_points[0].get_variable("x") == 0:
            mu0 = self.main_turbine.static_points[0].get_variable("visc")
            ss0 = self.main_turbine.static_points[0].get_variable("c")

        else:
            mu0, mu_l0, mu_g0 = self.evaluate_two_phase_visc(self.main_turbine.static_points[0].get_variable("p"),
                                                             self.main_turbine.static_points[0].get_variable("x"))
            ss0 = self.evaluate_two_phase_ss(self.main_turbine.static_points[0].get_variable("p"),
                                             self.main_turbine.static_points[0].get_variable("x"))

        Ma0 = first_speed / ss0
        h0 = self.main_turbine.static_points[0].get_variable("h")

        # Calculating Velocity and Thermodynamic Static Conditions Inside the Stator
        for i in range(self.options.n_stator):

            if i == 0:

                self.h[i] = self.main_turbine.static_points[0].get_variable("h")
                self.s[i] = self.main_turbine.static_points[0].get_variable("s")
                self.rho[i] = self.main_turbine.static_points[0].get_variable("rho")
                self.x[i] = self.main_turbine.static_points[0].get_variable("x")
                self.p[i] = self.main_turbine.static_points[0].get_variable("p")
                self.Ma[i] = Ma0
                self.mu[i] = mu0
                self.ss[i] = ss0
                self.epsilon[i], self.rho_l[i], self.rho_g[i] = self.evaluate_epsilon(self.p[i], self.x[i])
                self.chi[i], self.chi_s[i], self.mu_l[i], self.mu_g[i], self.mu[i] = self.evaluate_lockhart_martinelli(self.p[i], self.x[i])
                self.__stator_points[0] = self.main_turbine.static_points[0].duplicate()

            else:

                self.rho[i] = self.__stator_points[i-1].get_variable("rho")

                # Calculating Pressure Losses --> Churchill Correlation requires to know the Reynolds number, i.e. the
                # Lockhart-Martinelli parameter is used for retrieving the reynolds number and estimating the dpdl

                self.epsilon[i], self.rho_l[i], self.rho_g[i] = self.evaluate_epsilon(self.p[i-1], self.x[i-1])
                self.chi[i], self.chi_s[i], self.mu_l[i], self.mu_g[i], self.mu[i] = (
                                                           self.evaluate_lockhart_martinelli(self.p[i-1], self.x[i-1]))

                self.cosy[i] = 20

                self.phi_l[i] = 1 + (self.cosy[i] / self.chi[i]) + 1 / self.chi_s[i]
                self.phi_g[i] = self.chi[i] * self.phi_l[i]

                self.c1[i] = (1 - self.epsilon[i]) ** 4 * self.phi_l[i] ** 0.98
                self.c2[i] = self.epsilon[i] ** 4 * self.phi_g[i] ** 0.38

                self.dh_l[i] = np.sqrt(4 * A_eff[i] * (1 - self.epsilon[i]) / (self.c1[i] * np.pi))
                self.dh_g[i] = np.sqrt(4 * A_eff[i] * self.epsilon[i] / (self.c2[i] * np.pi))

                self.Re_l[i] = (self.m_dot_s * (1 - self.x[i - 1]) * self.dh_l[i] /
                                (A_eff[i] * (1 - self.epsilon[i]) * self.mu_l[i]))
                self.Re_g[i] = (self.m_dot_s * self.x[i - 1] * self.dh_g[i] /
                                (A_eff[i] * self.epsilon[i] * self.mu_g[i]))

                self.f_l[i] = self.churchill_friction_factor(self.Re_l[i], self.options.roughness, self.dh_l[i])
                self.dpdl_l[i] = (2 * self.f_l[i] * self.rho_l[i] *
                            (self.m_dot_s * self.x[i-1] / (self.rho_l[i] * A_eff[i] * (1 - self.epsilon[i]))) ** 2) / self.dh_l[i]
                self.dpdl_fr[i] = self.phi_l[i] * self.dpdl_l[i]

                self.p[i] = self.p[i-1] - dl[i] * self.dpdl_fr[i] - 0.5 * (self.rho[i] * (self.m_dot_s / (self.rho[i] * A_eff[i])) ** 2
                                                                           + self.rho[i-1] * (self.m_dot_s / (self.rho[i-1] * A_eff[i-1])) ** 2)
                self.h[i] = self.h[0] - 0.5 * (self.m_dot_s / (self.rho[i] * A_eff[i])) ** 2

                self.__stator_points[i].set_variable("p", self.p[i])
                self.__stator_points[i].set_variable("h", self.h[i])
                self.x[i] = self.__stator_points[i].get_variable("x")
