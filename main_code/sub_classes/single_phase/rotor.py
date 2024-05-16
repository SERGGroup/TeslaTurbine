import numpy as np

from main_code.base_classes.base_rotor import BaseRotor, BaseRotorStep, Speed, Position


class SPRotorStep(BaseRotorStep):

    def __init__(self, main_rotor, speed: Speed):

        super().__init__(main_rotor, speed)

    def get_variations(self, dr):

        rho = self.static_point.get_variable("rho")
        r_new = self.new_pos.r
        u = self.speed.u
        wt = self.speed.wt
        omega = self.new_pos.omega

        vr_new = self.m_dot /(2 * np.pi * self.geometry.b_channel * r_new * rho)

        dvr = vr_new - self.speed.vr

        coef = 6.5
        mu = self.static_point.get_variable("visc")
        ni = mu / rho

        wtFG = wt - dr * (-(10 / coef) * omega - (ni * 60 / (-vr_new * coef * self.geometry.b_channel ** 2) + 1 / r_new) * wt)
        wt_c = (wt + wtFG) / 2
        dwt = -dr *(-(10 / coef) * omega - (ni * 60/(-vr_new * coef * self.geometry.b_channel ** 2)+1 / r_new) * wt_c)
        wt_new = self.speed.wt + dwt

        dp = - dr * rho * ( omega ** 2 * r_new + 2 * wt_new * omega * coef / 6 + coef ** 2 / 30 * wt_new ** 2 / r_new - coef ** 2 / 30 *
                    vr_new * dvr / dr + 2 * coef * ni * vr_new / self.geometry.b_channel ** 2)

        dvt = dwt + self.main_rotor.omega * dr
        return dvr, dvt, dp, 0.

    def solve(self):

        pass


class SPRotor(BaseRotor):

    def __init__(self, main_turbine):

        super().__init__(main_turbine, SPRotorStep)

        self.isentropic_inlet = self.main_turbine.stator.isentropic_output
        self.intermediate_gap_point = self.main_turbine.points[0].duplicate()

    def evaluate_gap_losses(self):

        if self.gap_losses_control:

            # Evaluating Outlet Stator Pressure Loss
            ratio = self.main_turbine.geometry.stator.a_ratio_discharge

            ke = (1 - ratio) ** 2
            dh = self.main_turbine.points[0].get_variable("h") - self.isentropic_inlet.get_variable("h")
            rho1 = self.main_turbine.static_points[1].get_variable("rho")
            v_1ss = np.sqrt(2 * dh)
            v_1s = self.main_turbine.stator.phi_n * v_1ss
            DP_sbocco = 0.5 * ke * rho1 * v_1s ** 2

            # Evaluating Thermodynamic Conditions After Stator Loss
            P_01_R_star = self.main_turbine.points[1].get_variable("P") - DP_sbocco
            h_01_R_star = self.main_turbine.points[1].get_variable("h")               # Total Enthalpy Conservation
            h_1_R_star = self.main_turbine.static_points[1].get_variable("h")         # The Loss is Considered Iso-Enthalpic

            self.intermediate_gap_point.set_variable("h", h_01_R_star)
            self.intermediate_gap_point.set_variable("P", P_01_R_star)
            rho_1_R_star = self.intermediate_gap_point.get_variable("rho")

            # Evaluating Speed after Stator Loss
            vr_1_R_star = (self.main_turbine.stator.m_dot_s / self.geometry.n_channels) / (
                2 * np.pi * self.geometry.b_channel * ((self.geometry.d_out + 2 * self.geometry.gap) / 2) * rho_1_R_star)
            wr_1_R_star = vr_1_R_star

            # Evaluating Inlet Rotor Pressure Loss
            A_in_im = 2 * np.pi * (self.geometry.d_out / 2) * self.main_turbine.geometry.H_s
            A_out_im = self.geometry.n_channels * 2 * np.pi * (self.geometry.d_out / 2) * self.geometry.b_channel
            A2onA1 = (A_out_im / A_in_im)
            kc_1 = -0.12 * A2onA1 ** 4 + 1.02 * A2onA1 ** 3 - 1.28 * A2onA1 ** 2 - 0.12 * A2onA1 + 0.5
            DP_imbocco = 0.5 * kc_1 * rho_1_R_star * wr_1_R_star ** 2

            # Evaluating Thermodynamic Conditions After Rotor Loss
            P_01_R = P_01_R_star - DP_imbocco
            h_01_R = h_01_R_star
            h_1_R = h_1_R_star

            self.main_turbine.points[2].set_variable("P", P_01_R)
            self.main_turbine.points[2].set_variable("h", h_01_R)

            self.main_turbine.static_points[2].set_variable("P", P_01_R_star)
            self.main_turbine.static_points[2].set_variable("h", h_1_R)

            # Evaluating Speed after Rotor Loss
            vr_1_R = (self.main_turbine.stator.m_dot_s / self.geometry.n_channels) / (2 * np.pi * self.geometry.b_channel * (
                self.geometry.d_out / 2) * rho_1_R_star)
            v_1_R = np.sqrt(2 * (h_01_R - h_1_R))
            self.rotor_inlet_speed.init_from_codes("v", v_1_R, "v_r", vr_1_R)

            return self.rotor_inlet_speed

        else:
            self.main_turbine.points[1].copy_state_to(self.main_turbine.points[2])
            self.main_turbine.static_points[1].copy_state_to(self.main_turbine.static_points[2])
            self.rotor_inlet_speed = self.main_turbine.stator.speed_out
