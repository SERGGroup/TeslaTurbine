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

    d_ratio = 1.5
    throat_width = 0.0002
    Z_stat = 4
    H_s = 0.0005
    alpha1 = 85
    N_s = 50

    def __init__(self, main_geom):

        self.main_class = main_geom

    @property
    def d_int(self):
        return self.main_class.d_main

    @d_int.setter
    def d_int(self, d_int_in):
        self.main_class.d_main = d_int_in

    @property
    def r_int(self):
        return self.d_int / 2

    @property
    def chord(self):
        return 0.75 * self.r_int

    @property
    def alpha_rad(self):

        return self.alpha1 * np.pi / 180

class TeslaGeometry:

    d_main = 0.2

    def __init__(self):

        self.stator = StatorGeometry(self)
        self.rotor = RotorGeometry(self)
