import numpy as np


class BearingLoss:

    """
    This class is meant to be used for calculating the bearing moment losses.
    In the following, the assumption of SPHERICAL ROLLER BEARINGS has been made. Therefore, the parameters for the
    correlations are tuned for that specific kind of bearings. Further information can be found in the documentation directory
    """

    # Input Variables
    n = 0.
    nu = 0.
    d = 0.
    D = 0.
    dm = (d + D) / 2
    d_turb = 0.
    l_turb = 0.
    m_dot = 0.
    v = 0.
    p = 0.
    A_out = 0.


    # Spherical Roller Bearings Coefficients
    R1 = 1.6 * 10 ** (-6)
    R2 = 5.84
    R3 = 2.81 * 10 ** (-6)
    R4 = 5.8

    S1 = 3.62 * 10 ** (-3)
    S2 = 508
    S3 = 8.8 * 10 ** (-3)
    S4 = 117

    def __init__(self):

        self.M_rr = 0.
        self.M_sl = 0.
        self.M_seal = 0.
        self.M_drag = 0.
        self.M_lost = 0.
        self.n_bearings = 0

        # Radial Force is calculated by scaling the one of the lab-scale Tesla prototype
        self.F_r_ref = 2 * 9.81  # [N]
        self.D_ref = 0.12        # [m]
        self.L_ref = 0.1         # [m]

        self.F_r = (self.F_r_ref * (self.d_turb ** 2 / self.D_ref ** 2) * (self.l_turb / self.L_ref)) / self.n_bearings
        self.F_a = 0.5 * (self.m_dot * self.v + self.p * self.A_out)

    def __evaluate_mrr(self):

        """

        This method is used to evaluate the Rolling Frictional Resistant Torque of the bearing.
        It needs as input:
            - The rotating speed "n" (in rpm)
            - The internal / external diameter of the bearing "d" and "D" (in mm)
            - The kinematic viscosity of the lubricating oil "nu" (in mm2/s or cst)
            - Some empirical parameters depending on the type of bearings "Kz" and "Krs" (non-dimensional)

        """

        K_rs = 6 * 10 ** (-8)
        K_z = 0.

        phi_ihs = 1 / (1 + 1.84 * 10 ** (-9) * (self.n * self.dm) ** 1.28 * self.nu ** 0.64)
        phi_rs = 1 / np.exp(K_rs * self.nu * (self.d + self.D) * np.sqrt(K_z / (2 * (self.D - self.d))))

        G_rr = np.min(self.R1 * self.dm ** 1.85 * (self.F_r + self.R2 * self.F_a) ** 0.54, self.R3 * self.dm ** 2.3 * (
                self.F_r + self.R4 * self.F_a) ** 0.31)

        self.M_rr = phi_ihs * phi_rs * G_rr * (self.nu * self.n) ** 0.6

    def __evaluate_msl(self):

        """

        This method is used to evaluate the Sliding Frictional Resistant Torque of the bearing.
        It needs as input:
            - The rotating speed "n" (in rpm)
            - The internal / external diameter of the bearing "d" and "D" (in mm)
            - The kinematic viscosity of the lubricating oil "nu" (in mm2/s or cst)
            - Some empirical parameters depending on the type of bearings "mu_bl" and "mu_ehl" (non-dimensional)

        """
        mu_bl = 0.15
        mu_ehl = 0.1

        phi_bl = 1 / np.exp(2.6 * 10 ** (-8) * (self.nu * self.n) ** 1.4 * self.dm)
        mu_sl = phi_bl * mu_bl + (1 - phi_bl) * mu_ehl

        G_sl = np.min(self.S1 * self.dm ** 0.25 * (self.F_r ** 4 + self.S2 * self.F_a ** 4) ** 0.33, self.S3 * self.dm ** 0.94 * (self.F_r ** 3 + self.S4 * self.F_a ** 3) ** 0.33)

        self.M_sl = mu_sl + G_sl

    def __evaluate_mseal(self):

        pass

    def __evaluate_mdrag(self):

        pass

    def evaluate_bearing_losses(self):

        self.__evaluate_msl()
        self.__evaluate_mrr()
        self.__evaluate_mdrag()
        self.__evaluate_mseal()

        self.M_lost = self.M_rr + self.M_sl + self.M_seal + self.M_drag