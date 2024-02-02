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

    def theta_rel(self, theta_0=0., unit="rad"):

        if unit == "deg" or unit == "°":

            return self.theta * 180 / np.pi + theta_0

        return self.theta + theta_0

    def gamma_rel(self, gamma_0=0., unit="rad"):

        if unit == "deg" or unit == "°":
            return self.gamma * 180 / np.pi + gamma_0

        return self.gamma + gamma_0
