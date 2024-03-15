import numpy as np


class RotorGeometry:

    d_ratio = 3
    stator_gap = 0.0001
    b_channel = 0.0001
    roughness = 0.000001
    n_channels = 50

    def __init__(self, main_geom):

        self.main_class = main_geom

    @property
    def d_out(self):
        return self.main_class.d_main - 2 * self.stator_gap

    @d_out.setter
    def d_out(self, d_out_in):
        self.main_class.d_main = d_out_in + 2 * self.stator_gap

    @property
    def r_out(self):
        return self.d_out / 2

    @property
    def d_int(self):

        return self.d_out / self.d_ratio

    @property
    def r_int(self):
        return self.d_out / self.d_ratio / 2

    @property
    def dr_tot(self):

        return self.r_out - self.r_int


class StatorGeometry:

    d_ratio = 1.25
    Z_stat = 4
    H_s = 0.6
    N_s = 100

    # Those values need to be specified ONLY if working with a 0D Stator version
    throat_width = 0.003
    alpha1 = 85

    # Parameters for Stator Profiling
    t_u = 0.12      # [mm]
    ch_u = 12.5     # [mm]
    c_u = 0.19      # [mm]
    b_u = 0.17      # [mm]
    a_u = 0.1       # [mm]

    t_b = 0.12      # [mm]
    ch_b = 12.5     # [mm]
    c_b = 0.28      # [mm]
    b_b = 0.28      # [mm]
    a_b = 0.1       # [mm]

    u_0 = 12.5
    b_0 = 0.
    m_0 = 0.
    midline_0 = 6.25
    m_perp_0 = - 99999

    def __init__(self, main_geom):

        self.main_class = main_geom

    def get_prof(self, n, r_0, r_1, ax):

        ds = (r_0 - r_1) / (n-1)

        prof_array = np.empty((n, 5))
        prof_array[:] = np.nan

        for i in range(n):
            ax[i] = r_0 - (r_0 - i * ds)
            dup, dbottom, dmidline, m, m_perp = self.get_variation(ax[i])

            if i == 0:
                prof_array[i, 0] = self.u_0
                prof_array[i, 1] = self.b_0
                prof_array[i, 2] = self.midline_0
                prof_array[i, 3] = self.m_0
                prof_array[i, 4] = self.m_perp_0
            else:
                prof_array[i, 0] = self.u_0 + dup
                prof_array[i, 1] = self.b_0 + dbottom
                prof_array[i, 2] = self.midline_0 + dmidline
                prof_array[i, 3] = m
                prof_array[i, 4] = m_perp

        return prof_array

    def get_coordinates(self, n, r_0, r_1, ax):

        prof_array = self.get_prof(n, r_0, r_1, ax)
        r = np.zeros(n)
        m = np.zeros(n)
        X_SS = np.zeros(n)
        Y_SS = np.zeros(n)
        Y_PS = np.zeros(n)
        X_PS = np.zeros(n)
        X_LM = np.zeros(n)
        Y_LM = np.zeros(n)

        theta_SS = np.zeros(n)
        theta_PS = np.zeros(n)
        theta_LM = np.zeros(n)

        for i in range(n):

            r[i] = r_0 - ax[i]
            m[i] = prof_array[i, 3]
            X_SS[i] = prof_array[i,0]
            Y_SS[i] = np.sqrt(r[i] ** 2 - prof_array[i,0] ** 2)
            Y_PS[i] = np.sqrt(r[i] ** 2 - prof_array[i,1] ** 2)
            X_PS[i] = np.sqrt(r[i] ** 2 - Y_PS[i] ** 2)
            X_LM[i] = (X_SS[i] + X_PS[i]) / 2
            Y_LM[i] = (Y_SS[i] + Y_PS[i]) / 2

            theta_SS[i] = np.arcsin(X_SS[i] / r[i]) * 180 / np.pi
            theta_PS[i] = np.arcsin(X_PS[i] / r[i]) * 180 / np.pi
            theta_LM[i] = np.arcsin(X_PS[i] / r[i]) * 180 / np.pi

        alpha_out = (np.arcsin(X_LM[n-1] / r[n-1]) + np.arctan(m[n-1])) * 180 / np.pi

        return X_SS, Y_SS, Y_PS, X_PS, m, alpha_out

    def get_variation(self, ax):

        dup = 5 / 3 * self.t_u * self.ch_u * (self.a_u * (ax / self.ch_u) ** 2 + self.b_u * (ax / self.ch_u) ** 3 + self.c_u * (ax / self.ch_u) ** 4)
        dbottom = 5 / 3 * self.t_b * self.ch_b * (self.a_b * (ax / self.ch_b) ** 2 + self.b_b * (ax/ self.ch_b) ** 3 + self.c_b * (ax / self.ch_b) ** 4)
        dmidline = (dup + dbottom) / 2
        m =  5 / 6 * self.t_u * (2 * self.a_u * (ax / self.ch_u) + 3 * self.b_u * (ax / self.ch_u) ** 2 + 4 * self.c_u * (ax / self.ch_u) ** 3) + 5 / 6 * self.t_b * (2 * self.a_b * (ax / self.ch_b) + 3 * self.b_b * (ax/ self.ch_b) ** 2 + 4 * self.c_b * (ax / self.ch_b) ** 3)
        m_perp = - 1 / m

        return dup, dbottom, dmidline, m, m_perp

    def geometric_parameters(self, X_ss, Y_ss, Y_ps, X_ps, m, n):
        A_th = np.zeros(n)
        A_eff = np.zeros(n)

        for i in range(n):
            A_th[i] = np.sqrt((X_ss[i] - X_ps[i]) ** 2 + (Y_ss[i] - Y_ps[i]) ** 2) * self.H_s
            A_eff[i] = A_th[i] * np.cos(np.arctan(m[i]))

        throat_width = np.min(A_eff) / self.H_s

        return A_th, A_eff, throat_width

    @property
    def d_int(self):
        return self.main_class.d_main

    @d_int.setter
    def d_int(self, d_int_in):
        self.main_class.d_main = d_int_in

    @property
    def d_0(self):
        return self.main_class.d_main * self.d_ratio

    @property
    def r_0(self):
        return self.d_0 / 2

    @property
    def r_int(self):
        return self.d_int / 2

    @property
    def chord(self):
        return 0.75 * self.r_int

    @property
    def alpha_rad(self):

        return self.alpha1 * np.pi / 180


class BaseTeslaGeometry:

    d_main = 0.55

    def __init__(

            self,
            stator_geometry: type(StatorGeometry) = StatorGeometry,
            rotor_geometry: type(RotorGeometry) = RotorGeometry

    ):

        self.stator = stator_geometry(self)
        self.rotor = rotor_geometry(self)
