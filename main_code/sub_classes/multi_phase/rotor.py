import numpy as np
from main_code.base_classes.base_rotor import BaseRotor, BaseRotorStep, Speed
from .handlers.flow_pressure_losses_handler import flow_losses_handler
from .handlers.void_fraction_handler import void_fraction_handler

@flow_losses_handler
@void_fraction_handler
class TPRotorStep(BaseRotorStep):

    def __init__(self, main_rotor, speed: Speed):

        super().__init__(main_rotor, speed)

        self.liq_speed = Speed(speed.pos)
        self.vap_speed = Speed(speed.pos)

        self.liq_phase = self.main_rotor.input_point.duplicate()
        self.vap_phase = self.main_rotor.input_point.duplicate()

        self.chi_s = 0.
        self.chi = 0.
        self.Re_liq = 0.
        self.Re_gas = 0.
        self.dpdl_liq = 0.
        self.phi_liq = 0.
        self.D_liq = 0.

    def get_variations(self, dr):

        self.tp_evaluate_condition()
        self.evaluate_liq_pressure_losses()

        rho = self.thermo_point.get_variable("rho")
        r_new = self.new_pos.r
        vt = self.speed.vt
        vr = self.speed.vr
        beta = self.speed.beta

        vr_new = self.m_dot /(2 * np.pi * self.geometry.b_channel * r_new * rho)

        #
        # vtFG = vt - dr * (vt / r_new - self.phi_liq * self.dpdl_liq / (rho * vr_new) * np.sin(self.speed.beta))
        # vt_c = (vt + vtFG) / 2
        #

        dvr = vr_new - self.speed.vr
        dvt = dr * (self.speed.vt / r_new - self.phi_liq * self.dpdl_liq / (rho * vr_new) * np.sin(self.speed.beta))
        dP = - dr * ((rho * (vr ** 2 + vt ** 2) / r_new) - (self.phi_liq * self.dpdl_liq *
                np.cos(self.speed.beta)))

        return dvr, dvt, dP, 0.

    def tp_evaluate_condition(self):

        self.liq_phase.set_variable("P", self.thermo_point.get_variable("P"))
        self.vap_phase.set_variable("P", self.thermo_point.get_variable("P"))

        self.liq_phase.set_variable("x", 0)
        self.vap_phase.set_variable("x", 1)

        x = self.thermo_point.get_variable("x")
        m_g = self.m_dot * x
        m_l = self.m_dot * (1 - x)

        vr_l = m_l / (self.liq_phase.get_variable("rho") * (1 - self.epsilon) * 2 * np.pi * self.geometry.b_channel * self.pos.r)
        vr_g = m_g / (self.vap_phase.get_variable("rho") * self.epsilon * 2 * np.pi * self.geometry.b_channel * self.pos.r)

        self.liq_speed.init_from_codes("vr", vr_l, "beta", self.speed.beta)
        self.vap_speed.init_from_codes("vr", vr_g, "beta", self.speed.beta)

        self.chi_s = ((((1 - x) / x) ** 1.8) * ((self.vap_phase.get_variable("rho")) / (self.liq_phase.get_variable("rho")))
                 * ((self.vap_phase.get_variable("mu")) / (self.liq_phase.get_variable("mu"))) ** (-0.02))
        self.chi = np.sqrt(self.chi_s)

        self.Re_liq, self.Re_gas, self.D_liq, self.phi_liq = self.__iterate_cosy()

    def __iterate_cosy(self):

        cosy = 20
        n_it = 10

        # Initialization
        Re_liq = 0.
        Re_gas = 0.
        D_liq = 0.
        phi_liq = 0.

        for i in range(n_it):

            phi_liq = 1 + (cosy / self.chi ) + (1 / self.chi_s)
            phi_gas = phi_liq * self.chi

            c_1 = (1 - self.epsilon) ** 4 * phi_liq ** 0.98
            c_2 = self.epsilon ** 4 * phi_gas ** 0.38

            D_liq = np.sqrt(8 * self.pos.r * self.geometry.b_channel * (1 - self.epsilon) / c_1)
            D_gas = np.sqrt(8 * self.pos.r * self.geometry.b_channel * self.epsilon / c_2)

            Re_liq = (self.liq_phase.get_variable("rho") * self.liq_speed.w * D_liq) / self.liq_phase.get_variable("mu")
            Re_gas = (self.vap_phase.get_variable("rho") * self.vap_speed.w * D_gas) / self.vap_phase.get_variable("mu")

            if Re_liq > 4000 and Re_gas > 4000:
                cosyy = 20
            elif Re_liq < 4000 and Re_gas > 4000:
                cosyy = 12
            elif Re_liq > 4000 and Re_gas < 4000:
                cosyy = 10
            else:
                cosyy = 5

            if cosyy == cosy:
                break
            else:
                cosy = cosyy

        return Re_liq, Re_gas, D_liq, phi_liq

    def evaluate_liq_pressure_losses(self):

        A = (2.457 * np.log(1 / ((7 / self.Re_liq) ** 0.9 + 0.27 * self.geometry.roughness / self.D_liq))) ** 16
        B = (37530 / self.Re_liq) ** 16

        f = 2 * ((8 / self.Re_liq) ** 12 + 1 / ((A + B) ** 1.5)) ** (1 / 12)

        self.dpdl_liq = (2 * f * self.liq_phase.get_variable("rho") * self.liq_speed.w ** 2) / self.D_liq


@void_fraction_handler
class TPRotor(BaseRotor):

    def __init__(self, main_turbine):

        super().__init__(main_turbine, TPRotorStep)

        self.intermediate_point = self.main_turbine.points[0].duplicate()

        self.thermo_point = self.main_turbine.points[1].duplicate()

        self.liq_phase = self.main_turbine.points[1].duplicate()
        self.vap_phase = self.main_turbine.points[1].duplicate()

    def evaluate_gap_losses(self):

        """ The pressure losses in the stator-gap enlargement and in the gap-rotor contraction are considered to be
            iso-enthalpic --> both the total enthalpy and the static one are constant in the process, thus the
            absolute velocity is kept constant. Obviously, the pressure varies according to different correlations

            - ENLARGEMENT - J.G. Collier, J.R. Thome, Convective Boiling and Condensation, Clarendon Press, Oxford, 1994
            - CONTRACTION -  G.E. Geiger, Sudden Contraction Losses in Single and Two-Phase Flow

        """

        self.liq_phase.set_variable("P", self.main_turbine.static_points[1].get_variable("P"))
        self.vap_phase.set_variable("P", self.main_turbine.static_points[1].get_variable("P"))

        self.liq_phase.set_variable("x", 0)
        self.vap_phase.set_variable("x", 1)

        x = self.main_turbine.static_points[1].get_variable("x")
        # epsilon = 1 / (1 + (1 - x) / x * (self.vap_phase.get_variable("rho") / self.liq_phase.get_variable("rho")) ** (2 / 3))

        __tmp_point = list()

        for i in range(3):
            __tmp_point.append(self.main_turbine.points[0].duplicate())

        __tmp_point_liq1 = self.main_turbine.points[0].duplicate()
        __tmp_point_vap1 = self.main_turbine.points[0].duplicate()
        __tmp_point_liq2 = self.main_turbine.points[0].duplicate()
        __tmp_point_vap2 = self.main_turbine.points[0].duplicate()

        if self.gap_losses_control:

            A_in_sb = self.main_turbine.geometry.throat_width * self.main_turbine.geometry.H_s
            A_out_sb = (self.main_turbine.geometry.throat_width + self.geometry.gap / np.cos(
                np.radians(90 - self.main_turbine.geometry.alpha1))) * self.main_turbine.geometry.H_s
            A_ratio1 = A_in_sb / A_out_sb

            ke = (1 - A_in_sb / A_out_sb) ** 2
            rho1 = self.main_turbine.static_points[1].get_variable("rho")
            v_1s = self.main_turbine.stator.speed_out.v
            DP_sbocco = 0.5 * ke * rho1 * v_1s ** 2

            P_0 = self.main_turbine.points[1].get_variable("p") - DP_sbocco
            h_0 = self.main_turbine.points[1].get_variable("h")

            __tmp_point[1].set_variable("P", P_0)
            __tmp_point[1].set_variable("h", h_0)

            self.intermediate_point = __tmp_point[1].get_static_point(self.main_turbine.stator.speed_out.v)

            # Evaluation of Contraction Pressure Losses
            A_in_imb = 2 * np.pi * self.geometry.r_out * (self.main_turbine.geometry.disc_thickness + self.geometry.b_channel)
            A_out_imb = 2 * np.pi * self.geometry.r_out * self.geometry.b_channel
            A_ratio2 = A_out_imb / A_in_imb
            Cc = 1 - (1 - A_ratio2) / (2.08 * (1 - A_ratio2) + 0.5371)

            x_2 = self.intermediate_point.get_variable("x")

            n_it = 10

            self.main_turbine.points[2].set_variable("h", self.main_turbine.points[1].get_variable("h"))

            for n in range(n_it):

                __tmp_point[2].set_variable("x", x_2)
                __tmp_point[2].set_variable("h", self.main_turbine.static_points[2].get_variable("h"))

                P_g_2 = __tmp_point[2].get_variable("P")

                __tmp_point_liq2.set_variable("P", P_g_2)
                __tmp_point_vap2.set_variable("P", P_g_2)
                __tmp_point_vap2.set_variable("x", 1)
                __tmp_point_liq2.set_variable("x", 0)
                rho_l2 = __tmp_point_liq2.get_variable("rho")
                rho_v2 = __tmp_point_vap2.get_variable("rho")

                epsilon_star2 = 1 / (1 + (1 - x_2) / x_2 * (rho_v2 / rho_l2) ** (2 / 3))

                rho_first_guess2 = ((1 - x_2 ** 2) / (rho_l2 * (1 - epsilon_star2)) + x_2 ** 2 / (
                            rho_v2 * epsilon_star2)) ** (-1)
                rho_h_star2 = (x_2 / rho_v2 + (1 - x_2) / rho_l2) ** (-1)
                rho_sec_guess2 = ((1 - x_2) ** 3 / (rho_l2 * (1 - epsilon_star2)) ** 2 + x_2 ** 3 / (
                            rho_v2 * epsilon_star2) ** 2) ** (-1 / 2)
                G_con = self.m_dot_ch / A_out_imb
                DeltaP_con = G_con ** 2 * (rho_h_star2 / 2 * ((1 / (Cc * rho_sec_guess2) ** 2) - (A_ratio2 /
                                                                rho_sec_guess2) ** 2) + (1 - Cc) / rho_first_guess2) # RIVEDERE

                self.main_turbine.points[2].set_variable("P", P_0 - DeltaP_con)
                self.main_turbine.points[2].set_variable("h", self.main_turbine.points[1].get_variable("h"))

                self.main_turbine.static_points[2] = self.main_turbine.points[2].get_static_point(self.main_turbine.stator.speed_out.v)

                x2 = self.main_turbine.static_points[2].get_variable("x")

                if np.abs(x2 - x_2) < 0.000001:
                    break
                else:
                    x_2 = x2

            vr_R = self.m_dot_ch / (2 * np.pi * self.geometry.b_channel * (self.geometry.r_out) * __tmp_point[1].get_variable("rho"))
            v_R = self.main_turbine.stator.speed_out.v
            self.rotor_inlet_speed.init_from_codes("v", v_R, "v_r", vr_R)

            return self.rotor_inlet_speed

        else:

            self.main_turbine.points[1].copy_state_to(self.main_turbine.points[2])
            # self.main_turbine.static_points[1].copy_state_to(self.main_turbine.static_points[2])
            self.main_turbine.static_points[2].set_variable("x", self.main_turbine.static_points[1].get_variable("x"))
            self.main_turbine.static_points[2].set_variable("p", self.main_turbine.static_points[1].get_variable("p"))
            self.rotor_inlet_speed = self.main_turbine.stator.speed_out


    def evaluate_gap_losses_old(self):

        """ The pressure losses in the stator-gap enlargement and in the gap-rotor contraction are considered to be
            iso-enthalpic --> both the total enthalpy and the static one are constant in the process, thus the
            absolute velocity is kept constant. Obviously, the pressure varies according to different correlations

            - ENLARGEMENT - J.G. Collier, J.R. Thome, Convective Boiling and Condensation, Clarendon Press, Oxford, 1994
            - CONTRACTION -  G.E. Geiger, Sudden Contraction Losses in Single and Two-Phase Flow

        """

        self.liq_phase.set_variable("P", self.main_turbine.static_points[1].get_variable("P"))
        self.vap_phase.set_variable("P", self.main_turbine.static_points[1].get_variable("P"))

        self.liq_phase.set_variable("x", 0)
        self.vap_phase.set_variable("x", 1)

        __tmp_point = list()

        for i in range(3):
            __tmp_point.append(self.main_turbine.points[0].duplicate())

        __tmp_point_liq1 = self.main_turbine.points[0].duplicate()
        __tmp_point_vap1 = self.main_turbine.points[0].duplicate()
        __tmp_point_liq2 = self.main_turbine.points[0].duplicate()
        __tmp_point_vap2 = self.main_turbine.points[0].duplicate()

        if self.gap_losses_control:

            A_in_sb = self.main_turbine.geometry.throat_width * self.main_turbine.geometry.H_s
            # A_out_sb = ((self.main_turbine.geometry.throat_width / np.tan(
            #     np.radians(self.main_turbine.geometry.alpha1)) + self.geometry.gap / np.sin(
            #     np.radians(self.main_turbine.geometry.alpha1))) / np.cos(
            #     np.radians(self.main_turbine.geometry.alpha1)) - self.geometry.gap * np.tan(
            #     np.radians(self.main_turbine.geometry.alpha1)) - self.geometry.gap / np.tan(
            #     np.radians(self.main_turbine.geometry.alpha1))) * self.main_turbine.geometry.H_s
            A_out_sb = (self.main_turbine.geometry.throat_width + self.geometry.gap / np.cos(
                np.radians(90 - self.main_turbine.geometry.alpha1))) * self.main_turbine.geometry.H_s
            A_ratio1 = A_in_sb / A_out_sb

            # Evaluation of Enlargement Pressure Losses
            h_0_star = self.main_turbine.points[1].get_variable("h")
            h_star = self.main_turbine.static_points[1].get_variable("h")
            x_guess1 = self.main_turbine.static_points[1].get_variable("x")
            epsilon = 1 / (1 + (1 - x_guess1) / x_guess1 * (self.vap_phase.get_variable("rho") / self.liq_phase.get_variable("rho")) ** (2 / 3))

            rho_guess = ((1 - x_guess1 ** 2) / (self.liq_phase.get_variable("rho") * (1 - epsilon)) +
                         x_guess1 ** 2 / (self.vap_phase.get_variable("rho") * epsilon)) ** (-1)

            n_it = 500
            x_star = 0.
            P_0_star = 0.

            for n in range(n_it):

                __tmp_point[0].set_variable("x", x_guess1)
                __tmp_point[0].set_variable("h", h_star)

                P_g_star = __tmp_point[0].get_variable("P")

                __tmp_point_liq1.set_variable("P", P_g_star)
                __tmp_point_vap1.set_variable("P", P_g_star)
                __tmp_point_vap1.set_variable("x", 1)
                __tmp_point_liq1.set_variable("x", 0)
                rho_l = __tmp_point_liq1.get_variable("rho")
                rho_v = __tmp_point_vap1.get_variable("rho")

                epsilon_star = 1 / (1 + (1 - x_guess1) / x_guess1 * (rho_v / rho_l) ** (2 / 3))

                rho_first_guess = ((1 - x_guess1 ** 2) / (rho_l * (1 - epsilon)) + x_guess1 ** 2 / (rho_v * epsilon)) ** (-1)
                rho_h_star = (x_guess1 / rho_v + (1 - x_guess1) / rho_l) ** (-1)
                rho_sec_guess = ((1 - x_guess1) ** 3 / (rho_l * (1 - epsilon)) ** 2 + x_guess1 ** 3 / (rho_v * epsilon) ** 2) ** (-1/2)
                G_en = self.main_turbine.stator.m_dot_s / (A_in_sb * self.main_turbine.stator.geometry.Z_stat)
                DeltaP_en = G_en ** 2 / (2 * rho_l) * (2 * rho_l * A_ratio1 * (A_ratio1 - 1) / rho_first_guess - rho_h_star * rho_l * (A_ratio1 - 1) / rho_sec_guess ** 2)
                P_0_star = self.main_turbine.points[1].get_variable("p") + DeltaP_en

                __tmp_point[1].set_variable("P", P_0_star)
                __tmp_point[1].set_variable("h", h_0_star)

                # Overwriting Thermodynamic Intermediate Point
                self.intermediate_point.set_variable("h", h_0_star - 0.5 * self.main_turbine.stator.speed_out.v ** 2)
                self.intermediate_point.set_variable("p", P_0_star - 0.5 * rho_first_guess * self.main_turbine.stator.speed_out.v ** 2)

                x_star = self.intermediate_point.get_variable("x")
                p_star = self.intermediate_point.get_variable("p")

                __tmp_point_liq1.set_variable("P", p_star)
                __tmp_point_vap1.set_variable("P", p_star)
                __tmp_point_vap1.set_variable("x", 1)
                __tmp_point_liq1.set_variable("x", 0)
                rho_l1 = __tmp_point_liq1.get_variable("rho")
                rho_v1 = __tmp_point_vap1.get_variable("rho")

                epsilon_out = 1 / (1 + (1 - x_star) / x_star * (rho_v1 / rho_l1) ** (2 / 3))

                # if abs(epsilon_out - epsilon) < 0.00001:
                #     break
                # elif epsilon_out < epsilon:
                #     x_guess1 = x_guess1 + (x_star + x_guess1) / 100
                # else:
                #     x_guess1 = x_guess1 - (x_star + x_guess1) / 100

                if np.abs(x_star - x_guess1) < 0.001:
                    break
                elif x_star > x_guess1:
                    x_guess1 = x_guess1 + (x_star + x_guess1) / 500
                else:
                    x_guess1 = x_guess1 - (x_star + x_guess1) / 500


                # else:
                #      x_guess1 = x_star


                n_out = n + 1

            P_int = P_0_star
            # Evaluation of Contraction Pressure Losses
            A_in_imb = 2 * np.pi * self.geometry.r_out * (self.main_turbine.geometry.disc_thickness + self.geometry.b_channel)
            A_out_imb = 2 * np.pi * self.geometry.r_out * self.geometry.b_channel
            A_ratio2 = A_out_imb / A_in_imb
            Cc = 1 - (1 - A_ratio2) / (2.08 * (1 - A_ratio2) + 0.5371)

            x_guess2 = x_star
            n2 = 0

            self.main_turbine.points[2].set_variable("h", self.main_turbine.points[1].get_variable("h"))

            for n in range(n_it):

                __tmp_point[2].set_variable("x", x_guess2)
                __tmp_point[2].set_variable("h", self.main_turbine.static_points[2].get_variable("h"))

                P_g_2 = __tmp_point[2].get_variable("P")

                __tmp_point_liq2.set_variable("P", P_g_2)
                __tmp_point_vap2.set_variable("P", P_g_2)
                __tmp_point_vap2.set_variable("x", 1)
                __tmp_point_liq2.set_variable("x", 0)
                rho_l2 = __tmp_point_liq2.get_variable("rho")
                rho_v2 = __tmp_point_vap2.get_variable("rho")

                epsilon_star2 = 1 / (1 + (1 - x_guess2) / x_guess2 * (rho_v2 / rho_l2) ** (2 / 3))

                rho_first_guess2 = ((1 - x_guess2 ** 2) / (rho_l2 * (1 - epsilon_star2)) + x_guess2 ** 2 / (
                            rho_v2 * epsilon_star2)) ** (-1)
                rho_h_star2 = (x_guess2 / rho_v2 + (1 - x_guess2) / rho_l2) ** (-1)
                rho_sec_guess2 = ((1 - x_guess2) ** 3 / (rho_l2 * (1 - epsilon_star2)) ** 2 + x_guess2 ** 3 / (
                            rho_v2 * epsilon_star2) ** 2) ** (-1 / 2)
                G_con = self.m_dot_ch / A_out_imb
                DeltaP_con = G_con ** 2 * (rho_h_star2 / 2 * ((1 / (Cc * rho_sec_guess2) ** 2) - (A_ratio2 /
                                                                rho_sec_guess2) ** 2) + (1 - Cc) / rho_first_guess2) # RIVEDERE

                self.main_turbine.points[2].set_variable("P", P_int - DeltaP_con)
                self.main_turbine.points[2].set_variable("h", self.main_turbine.points[1].get_variable("h"))

                self.main_turbine.static_points[2].set_variable("h",
                                                                self.main_turbine.static_points[1].get_variable("h"))
                self.main_turbine.static_points[2].set_variable("s", self.main_turbine.points[2].get_variable("s"))

                x2 = self.main_turbine.static_points[2].get_variable("x")
                p2 = self.main_turbine.static_points[2].get_variable("P")
                if np.abs(x2 - x_guess2) < 0.000001:
                    break
                else:
                    x_guess2 = x2

            vr_R = self.m_dot_ch / (2 * np.pi * self.geometry.b_channel * (self.geometry.r_out) * __tmp_point[1].get_variable("rho"))
            v_R = self.main_turbine.stator.speed_out.v
            self.rotor_inlet_speed.init_from_codes("v", v_R, "v_r", vr_R)

            return self.rotor_inlet_speed

        else:

            self.main_turbine.points[1].copy_state_to(self.main_turbine.points[2])
            # self.main_turbine.static_points[1].copy_state_to(self.main_turbine.static_points[2])
            self.main_turbine.static_points[2].set_variable("x", self.main_turbine.static_points[1].get_variable("x"))
            self.main_turbine.static_points[2].set_variable("p", self.main_turbine.static_points[1].get_variable("p"))
            self.rotor_inlet_speed = self.main_turbine.stator.speed_out
