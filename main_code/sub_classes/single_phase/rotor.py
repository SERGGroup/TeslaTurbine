import numpy as np

from main_code.base_classes.base_rotor import BaseRotor, BaseRotorStep, Speed, Position


class SPRotorStep(BaseRotorStep):

    def __init__(self, main_rotor, speed: Speed):

        super().__init__(main_rotor, speed)

    def get_variations(self, dr):

        rho = self.thermo_point.get_variable("rho")
        r_new = self.new_pos.r
        u = self.speed.u
        wt = self.speed.wt
        omega = self.new_pos.omega

        # Ricontrollare questa portata m_dot!!
        vr_new = self.m_dot /(2 * np.pi * self.geometry.b_channel * r_new * rho)

        dvr = vr_new - self.speed.vr

        coef = 6.5
        mu = self.thermo_point.get_variable("visc")
        ni = mu / rho

        wtFG = wt - dr * (-(10 / coef) * omega - (ni * 60 / (-vr_new * coef * self.geometry.b_channel ** 2) + 1 / r_new) * wt)
        wt_c = (wt + wtFG) / 2
        dwt = -dr *(-(10 / coef) * omega - (ni * 60/(-vr_new * coef * self.geometry.b_channel ** 2)+1 / r_new) * wt_c)
        wt_new = self.speed.wt + dwt

        dP = - dr * rho * ( omega ** 2 * r_new + 2 * wt_new * omega * coef / 6 + coef ** 2 / 30 * wt_new ** 2 / r_new - coef ** 2 / 30 *
                    vr_new * dvr / dr + 2 * coef * ni * vr_new / self.geometry.b_channel ** 2)

        return dvr, dwt, dP, 0.

    def solve(self):

        pass


class SPRotor(BaseRotor):

    def __init__(self, main_turbine):

        super().__init__(main_turbine, SPRotorStep)

        self.isentropic_inlet = self.main_turbine.stator.isentropic_output
        self.intermediate_gap_point = self.main_turbine.points[0].duplicate()

        inlet_pos = Position(self.geometry.r_out, 0)
        self.rotor_inlet_speed = Speed(inlet_pos)

    def solve(self):

        self.rotor_points = list()
        self.rotor_inlet_speed = self.evaluate_gap_losses()

        self.omega_in = self.rotor_inlet_speed.vt / ((self.dv_perc + 1) * self.geometry.r_out)

        first_pos = Position(self.geometry.r_out, self.omega_in)
        first_speed = Speed(position=first_pos)

        # TODO: Change to static pressure model
        first_speed.equal_absolute_speed_to(self.rotor_inlet_speed)

        self.rothalpy = (self.main_turbine.static_points[2].get_variable("h") + (first_speed.w ** 2) / 2 -
                         (first_speed.u ** 2) / 2)

        first_step = self.rotor_step_cls(self, first_speed)

        new_step = first_step

        #
        # For detailed explanation on rotor discretization check:
        #
        #   "main_code/base_classes/other/rotor discretization explaination.xlsx"
        #

        dr_tot = self.geometry.dr_tot
        b = self.options.integr_variable * dr_tot
        a = np.power(dr_tot / b + 1, 1 / self.options.n_rotor)

        for i in range(self.options.n_rotor):

            if self.options.profile_rotor:
                self.rotor_points.append(new_step)

            dr = a ** i * (a - 1) * b
            new_step = new_step.get_new_step(dr)

        if self.options.profile_rotor:
            self.rotor_points.append(new_step)

    def evaluate_gap_losses(self):

        # Evaluating Outlet Stator Pressure Loss
        A_in_sb = self.main_turbine.stator.geometry.throat_width * self.main_turbine.stator.geometry.H_s
        A_out_sb = ((self.main_turbine.stator.geometry.throat_width / np.tan(np.radians(90 - self.geometry.alpha1)) + self.geometry.gap / np.sin(
            np.radians(90 - self.geometry.alpha1))) / np.cos(np.radians(90 - self.geometry.alpha1)) -
                    self.geometry.gap * np.tan(np.radians(90 - self.geometry.alpha1)) - self.geometry.gap / np.tan(
                    np.radians(90 - self.geometry.alpha_1PS))) * self.main_turbine.stator.geometry.H_s

        rapporto = self.geometry.Z_stat * A_out_sb / (2 * np.pi * (self.geometry.d_out / 2) * self.main_turbine.stator.geometry.H_s)
        ke = (1 - A_in_sb / A_out_sb) ** 2
        dh = self.main_turbine.points[0].get_variable("h") - self.isentropic_inlet.get_variable("h")
        rho1 = self.main_turbine.static_points[1].get_variable("rho")
        v_1s = np.sqrt(2 * dh)
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
        A_in_im = 2 * np.pi * (self.geometry.d_out / 2) * self.geometry.H_s
        A_out_im = self.geometry.n_discs * 2 * np.pi * (self.geometry.d_out / 2) * self.geometry.b_channel
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

