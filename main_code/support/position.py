import numpy as np


class Position:

    def __init__(self, r, omega, theta=0, gamma=0, t=0):

        self.t = t

        self.r = r
        self.theta = theta
        self.gamma = gamma

        self.__omega = omega

    @property
    def rpm(self):

        return self.__omega / np.pi * 30

    @property
    def omega(self):

        return self.__omega

    @property
    def u(self):

        return self.__omega * self.r
