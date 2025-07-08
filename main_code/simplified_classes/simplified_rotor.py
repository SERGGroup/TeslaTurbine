from main_code.base_classes.base_rotor import BaseRotorStep, BaseRotor
from main_code.base_classes.support import Speed
from abc import ABC, abstractmethod
import numpy as np

class SimplifiedBaseRotorStep(BaseRotorStep, ABC):

    def get_new_step(self, dr):

        # Update Position
        self.new_pos = self.speed.get_new_position(dr)

        dvr, dvt, dp, dh = self.get_variations(-dr)

        # Update Speed
        new_speed = Speed(self.new_pos)
        vt_new = self.speed.vt + dvt
        vr_new = self.speed.vr + dvr
        new_speed.init_from_codes("vt", vt_new, "vr", vr_new)

        # Init a new step
        new_class = self.self_class()
        new_step = new_class(main_rotor=self.main_rotor, speed=new_speed)

        # Update Thermodynamic Point
        p_new = self.thermo_point.get_variable("P") + dp
        h_new = self.thermo_point.get_variable("H") + dh

        new_step.thermo_point.set_variable("P", p_new)
        new_step.thermo_point.set_variable("H", h_new)

        return new_step

    def get_variations(self, dr):
        """ Function to evaluate the variations of the thermodynamic variables
        and speed given the discrimination of the Navier-Stokes equations.
        Derivatives of thermodynamic variables are retrieved from the REFPROP.
        The derivatives of rothalpy and the friction factor should be provided
        by the user
        """

        r = self.pos.r
        b_chnl = self.geometry.b_channel
        vr = self.speed.vr
        vt = self.speed.vt
        wt = self.speed.wt
        om = self.pos.omega

        drho_dp = self.thermo_point.get_derivative("rho", "P", "H")
        drho_dh = self.thermo_point.get_derivative("rho", "H", "P")
        rho = self.thermo_point.get_variable("rho")

        d_hyd = 2 * (b_chnl * 2 * np.pi * r) / (b_chnl + 2 * np.pi * r)
        f = abs((self.get_friction_coefficient() * self.speed.w ** 2) / (2 * d_hyd))
        di = self.get_heat_losses() / (b_chnl * rho * vr)

        dvt = f * self.speed.sin_beta / vr - vt / r

        A = di + om * vt - wt * dvt
        B = f * self.speed.cos_beta + vt ** 2 / r
        C = 1 - vr ** 2 * (drho_dp + drho_dh / rho)
        D = r * (rho * drho_dp * B + drho_dh * A) / rho

        dvr = - (1 + D) / C * vr / r
        dp = rho * (B - vr * dvr)
        dh = A - vr * dvr

        return dvr * dr, dvt * dr, dp * dr, dh * dr

    @abstractmethod
    def get_friction_coefficient(self):
        """

            This function must return the friction coefficient for the current
            rotor step. Darcy Factor should be returned.

        """
        return 0.

    @abstractmethod
    def get_heat_losses(self):
        """

            This function must return the heat losses for the current rotor step
            in [W/m^2]. It is used to evaluate the rothalpy balance.
            For adiabatic turbines, return 0.

        """
        return 0.


class SimplifiedRotorStep(SimplifiedBaseRotorStep):
    """
        This class is a simplified version of the rotor step, which is used to
        evaluate the performance of the rotor in a Tesla turbine. It is a
        simplified version of the BaseRotorStep class, which is used to evaluate
        the performance of the rotor in a Tesla turbine.
    """

    def get_friction_coefficient(self):

        """Calculate the friction coefficient based on the Churchill Equation"""
        r = self.pos.r
        b_chnl = self.geometry.b_channel
        d_hyd = 2 * (b_chnl * 2 * np.pi * r) / (b_chnl + 2 * np.pi * r)
        rho = self.thermo_point.get_variable("rho")
        mu = self.thermo_point.get_variable("mu")
        re = rho * np.abs(self.speed.w) * d_hyd / mu

        A = (2.457 * np.log(1 / ((7 / re) ** 0.9 + 0.27 * self.geometry.roughness / d_hyd))) ** 16
        B = (37530 / re) ** 16
        return 2 * ((8 / re) ** 12 + 1 / ((A + B) ** 1.5)) ** (1 / 12)

    def get_heat_losses(self):
        return 0.


class SimplifiedRotor(BaseRotor):

    def __init__(self, main_turbine):

        super().__init__(main_turbine, SimplifiedRotorStep)

        self.isentropic_inlet = self.main_turbine.stator.isentropic_output
        self.intermediate_gap_point = self.main_turbine.points[0].duplicate()

    def evaluate_gap_losses(self):

        if self.gap_losses_control:

            # Evaluating Outlet Stator Pressure Loss
            A_in_sb = self.main_turbine.geometry.throat_width * self.main_turbine.geometry.H_s
            # A_out_sb = ((self.main_turbine.geometry.throat_width / np.tan(np.radians(90 - self.main_turbine.geometry.alpha1)) + self.geometry.gap / np.sin(
            #     np.radians(90 - self.main_turbine.geometry.alpha1))) / np.cos(np.radians(90 - self.main_turbine.geometry.alpha1)) -
            #         self.geometry.gap * np.tan(np.radians(90 - self.main_turbine.geometry.alpha1)) - self.geometry.gap / np.tan(
            #         np.radians(90 - self.geometry.alpha_1PS))) * self.main_turbine.geometry.H_s
            A_out_sb = (self.main_turbine.geometry.throat_width + self.geometry.gap / np.cos(np.radians(90 - self.main_turbine.geometry.alpha1))) * self.main_turbine.geometry.H_r

            ''' Borda-Carnot Loss Coefficient --- Basic Hypothesis
            - incompressible flow (or, at least, very low Mach number)
            - subsonic flow
            '''

            ke = (1 - A_in_sb / A_out_sb) ** 2

            ''' Nouri-Burujerdi Loss Coefficient --- Basic Hypothesis
            - compressible flow
            - Mach number < 0.6
            '''

            # M = 0.6
            # ke = 1.06 * (1 - A_in_sb / A_out_sb) ** 1.9 * np.exp((3 * M ** 2.87) / ((1 - A_in_sb / A_out_sb) ** 0.757))

            dh = self.main_turbine.points[0].get_variable("h") - self.isentropic_inlet.get_variable("h")
            rho1 = self.main_turbine.static_points[1].get_variable("rho")
            v_1ss = np.sqrt(2 * dh)
            v_1s = self.main_turbine.stator.phi_n * v_1ss
            DP_sbocco = 0.5 * ke * rho1 * v_1s ** 2

            # Evaluating Thermodynamic Conditions After Stator Loss
            P_01_R_star = self.main_turbine.points[1].get_variable("P") - DP_sbocco
            h_01_R_star = self.main_turbine.points[1].get_variable("h")               # Total Enthalpy Conservation
            h_1_R_star = h_01_R_star - 0.5 * v_1s ** 2         # The Loss is Considered Iso-Enthalpic

            self.intermediate_gap_point.set_variable("h", h_01_R_star)
            self.intermediate_gap_point.set_variable("P", P_01_R_star)
            rho_1_R_star = self.intermediate_gap_point.get_variable("rho")

            # Evaluating Speed after Stator Loss
            # vr_1_R_star = (self.main_turbine.stator.m_dot_s / self.geometry.n_channels) / (
            #     2 * np.pi * self.geometry.b_channel * ((self.geometry.d_out + 2 * self.geometry.gap) / 2) * rho_1_R_star)
            vr_1_R_star = (self.main_turbine.stator.m_dot_s / self.main_turbine.n_packs) / (
                2 * np.pi * self.geometry.b_channel * ((self.geometry.d_out + 2 * self.geometry.gap) / 2) * rho_1_R_star)
            wr_1_R_star = vr_1_R_star

            # Evaluating Inlet Rotor Pressure Loss
            A_in_im = 2 * np.pi * (self.geometry.d_out / 2) * self.main_turbine.geometry.H_r
            A_out_im = self.geometry.n_discs * 2 * np.pi * (self.geometry.d_out / 2) * self.geometry.b_channel
            A2onA1 = (A_out_im / A_in_im)
            kc_1 = -0.12 * A2onA1 ** 4 + 1.02 * A2onA1 ** 3 - 1.28 * A2onA1 ** 2 - 0.12 * A2onA1 + 0.5
            DP_imbocco = 0.5 * kc_1 * rho_1_R_star * wr_1_R_star ** 2

            # Evaluating Thermodynamic Conditions After Rotor Loss
            P_01_R = P_01_R_star - DP_imbocco
            h_01_R = h_01_R_star
            h_1_R = h_01_R - 0.5 * v_1s ** 2

            self.main_turbine.points[2].set_variable("P", P_01_R)
            self.main_turbine.points[2].set_variable("h", h_01_R)

            self.main_turbine.static_points[2].set_variable("P", P_01_R - 0.5 * rho_1_R_star * v_1s ** 2)
            self.main_turbine.static_points[2].set_variable("h", h_1_R)

            # Evaluating Speed after Rotor Loss
            # vr_1_R = (self.main_turbine.stator.m_dot_s / self.geometry.n_channels) / (2 * np.pi * self.geometry.b_channel * (
            #     self.geometry.d_out / 2) * rho_1_R_star)
            vr_1_R = (self.main_turbine.stator.m_dot_s / self.main_turbine.n_packs) / (2 * np.pi * self.geometry.b_channel * (
                self.geometry.d_out / 2) * rho_1_R_star)
            v_1_R = np.sqrt(2 * (h_01_R - h_1_R))
            self.rotor_inlet_speed.init_from_codes("v", v_1_R, "v_r", vr_1_R)

            return self.rotor_inlet_speed

        else:
            self.main_turbine.points[1].copy_state_to(self.main_turbine.points[2])
            self.main_turbine.static_points[1].copy_state_to(self.main_turbine.static_points[2])
            self.rotor_inlet_speed = self.main_turbine.stator.speed_out