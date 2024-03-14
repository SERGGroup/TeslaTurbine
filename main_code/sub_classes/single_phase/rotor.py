import numpy as np

from main_code.base_classes.base_rotor import BaseRotor, BaseRotorStep, Speed

class SPRotorStep(BaseRotorStep):

    def __init__(self, main_rotor, speed: Speed):

        super().__init__(main_rotor, speed)

    def get_variations(self, dr):

        rho = self.thermo_point.get_variable("rho")
        r_new = self.new_pos.r
        u = self.new_pos.u
        wt = self.new_pos.wt

        vr_new = self.m_dot /(2 * np.pi * self.geometry.b_channel * r_new * rho)

        dvr = vr_new - self.speed.vr

        coef = 6.5
        mu1 = self.static_output_point.get_variable("visc")
        n1 = mu1 / rho

        wtFG = wt - dr * (
                -(10 / coef) * self.__omega - (n1 * 60 / (-vr_new * coef * self.geometry.b_channel ^ 2) + 1 / r_new) * wt)
        wt_c = (wt + wtFG) / 2

        dwt = -dr *(-(10 / coef) * self.__omega - (n1 * 60/(-vr_new * coef * self.geometry.b_channel ^ 2)+1 / r_new) * wt_c)
        wt_new = self.speed.wt + dwt


        dP =  - dr * rho * ( self.__omega ^ 2 * r_new + 2 * wt_new * self.__omega * coef / 6 + coef ^ 2 / 30 * wt_new ^ 2 / r_new - coef ^ 2 / 30 *
                    vr_new * dvr / dr + 2 * coef * n1 * vr_new / self.geometry.b_channel ^ 2)

        return dvr / dr, dwt / dr, dP / dr, 0.

    def solve(self):

        pass

class SPRotor(BaseRotor):

    def __init__(self, main_turbine):

        super().__init__(main_turbine, SPRotorStep)

    def evaluate_gap_losses(self):
        A_in_sb = self.geometry.throat_width * self.geometry.H_s
        A_out_sb = ((self.geometry.throat_width / np.tan(np.radians(90 - self.geometry.alpha1)) + self.geometry.GAP / np.sin(
            np.radians(90 - self.geometry.alpha1))) / np.cos(np.radians(90 - self.geometry.alpha1)) -
                    self.geometry.GAP * np.tan(np.radians(90 - self.geometry.alpha1)) - self.geometry.GAP / np.tan(
                    np.radians(90 - self.geometry.alpha_1PS))) * self.geometry.H_s

        rapporto = self.geometry.Z_stat * A_out_sb / (2 * np.pi * self.geometry.r1_s * self.geometry.H_s)
        ke = (1 - A_in_sb / A_out_sb) ** 2
        dh = self.input_point.get_variable("h") - self.__tmp_points[0].get_variable("h")
        rho1 = self.static_output_point.get_variable("rho")
        v_1s = np.sqrt(2 * dh)

        DP_sbocco = 0.5 * ke * rho1 * v_1s ** 2
        P_01_R_star = self.p_out - DP_sbocco
        h_01_R_star = self.input_point.get_variable("h")
        h_1_R_star = self.static_output_point.set_variable("h", self.output_point.get_variable("h"))
        self.__tmp_points[0].set_variable("p", self.p_out)


        rho_1_R_star = self.static_output_point.get_variable("rho")
        vr_1_R_star = self.m_dot / (2 * np.pi * self.geometry.b * self.geometry.r1_s * rho_1_R_star)
        wr_1_R_star = vr_1_R_star


        A_in_im = 2 * np.pi * self.geometry.r1_s * self.geometry.H_s
        A_out_im = self.geometry.n_discs * 2 * np.pi * self.geometry.r1_s * self.geometry.b
        A2onA1 = (A_out_im / A_in_im)
        kc_1 = -0.12 * A2onA1 ** 4 + 1.02 * A2onA1 ** 3 - 1.28 * A2onA1 ** 2 - 0.12 * A2onA1 + 0.5
        DP_imbocco = 0.5 * kc_1 * rho_1_R_star * wr_1_R_star ** 2
        P_01_R = P_01_R_star - DP_imbocco
        h_01_R = h_01_R_star
        h_1_R = h_1_R_star



        self.main_turbine.points[2].set_variable("P", P_01_R)
        self.main_turbine.points[2].set_variable("h", h_01_R)

        self.main_turbine.static_points[2].set_variable("P", P_01_R_star)
        self.main_turbine.static_points[2].set_variable("h", h_1_R)
