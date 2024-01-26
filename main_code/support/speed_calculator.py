import numpy as np


class Position:

    def __init__(self, r, theta, gamma, omega):

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

class Speed:

    __possible_codes = ["v", "v_t", "v_r", "w", "w_t", "w_r", "alpha", "beta"]

    __v = 0.
    __vt = 0.
    __vr = None

    __w = None
    __wt = None
    __wr = None

    __alpha = None
    __beta = None

    def __init__(self, position: Position):

        self.pos = position

    def init_from_codes(self, first_code, first_value, second_code, second_value):

        """
        Input main_code can be:

            - v ................ absolute speed magnitude [m/s]
            - v_t .............. absolute tangential speed magnitude [m/s]
            - v_r .............. absolute radial speed magnitude [m/s]

            - w ................ relative speed magnitude [m/s]
            - w_t .............. relative tangential speed magnitude [m/s]
            - w_r .............. relative radial speed magnitude [m/s]

            - alpha ............ absolute speed angle (to radial) [m/s]
            - beta ............. relative speed angle (to radial) [m/s]

        To allow for the calculation to begin you should provide two independent parameters.
        The following table represent the dependencies between parameters:

        |param     |  v  |  v_t  |  v_r  |  w  |  w_t  |  w_r  |   beta  |  alpha  |
        --------------------------------------------------------------------------
        |v         |  -  |   v   |   v   |  v  |   v   |   v   |    -    |    v    |
        |v_t       |  v  |   -   |   v   |  v  |   -   |   v   |    v    |    v    |
        |v_r       |  v  |   v   |   -   |  v  |   v   |   -   |    v    |    v    |
        |w         |  v  |   v   |   v   |  -  |   v   |   v   |    v    |    -    |
        |w_t       |  v  |   -   |   v   |  v  |   -   |   v   |    v    |    v    |
        |w_r       |  v  |   v   |   -   |  v  |   v   |   -   |    v    |    v    |
        |beta      |  -  |   v   |   v   |  v  |   v   |   v   |    -    |    -    |
        |alpha     |  v  |   v   |   v   |  -  |   v   |   v   |    -    |    -    |

        meaning that the only dependent variables are:

            - alpha and beta
            - v_t and w_t
            - v_r and w_r

        in addition, other two couples of parameter cannot be specified as they are
        not enough to fully define the system:

            - alpha and w
            - beta and v

        :param second_value: Float, variable
        :param second_code: String Input Code
        :param first_value: Float, variable
        :param first_code: String Input Code
        :return: None
        """

        if self.__codes_are_independent(first_code, second_code):

            self.__evaluate(first_code, first_value, second_code, second_value)
            self.__evaluate_other()

    def get_new_position(self, dt):

        new_r = self.pos.r + self.__vr * dt
        new_theta = self.pos.theta + self.__vt / self.pos.r * dt
        new_gamma = self.pos.gamma + self.__wt / self.pos.r * dt

        return Position(new_r, new_theta, new_gamma, self.pos.omega)

    @property
    def v(self):

        return self.__v

    @property
    def vt(self):

        return self.__vt

    @property
    def vr(self):

        return self.__vr

    @property
    def w(self):

        return self.__w

    @property
    def wt(self):

        return self.__wt

    @property
    def wr(self):

        return self.__wr

    @property
    def alpha(self):

        return self.__alpha

    @property
    def beta(self):

        return self.__beta

    def __check_code_correct(self, code):

        return code in self.__possible_codes

    def __evaluate(self, first_code, first_value, second_code, second_value):

        input_codes = [first_code, second_code]
        input_values = [first_value, second_value]

        if "v" in input_codes:

            i = input_codes.index("v")
            self.__v = input_values[i]

            if "v_t" in input_codes or "w_t" in input_codes:

                if "w_t" in input_codes:

                    i = input_codes.index("w_t")
                    self.__wt = input_values[i]
                    self.__vt = self.__wt + self.pos.u

                else:

                    i = input_codes.index("v_t")
                    self.__vt = input_values[i]

            elif "v_r" in input_codes or "w_r" in input_codes:

                if "w_r" in input_codes:

                    i = input_codes.index("w_r")

                else:

                    i = input_codes.index("v_r")

                self.__vr = input_values[i]
                self.__wr = self.__vr
                self.__vt = np.sqrt(self.__v ** 2 - self.__vr ** 2)

            elif "w" in input_codes:

                i = input_codes.index("w")
                self.__w = input_values[i]

                u = self.pos.u
                self.__vt = 1 / 2 * ((self.__v ** 2 - self.__w ** 2) / u + u)

            elif "alpha" in input_codes:

                    i = input_codes.index("alpha")
                    self.__alpha = input_values[i]
                    self.__vt = self.__v * np.sin(self.__alpha)

        elif "v_t" in input_codes or "w_t" in input_codes:

            if "v_t" in input_codes:

                i = input_codes.index("v_t")
                self.__vt = input_values[i]
                self.__wt = self.__vt - self.pos.u

            else:

                i = input_codes.index("v_t")
                self.__wt = input_values[i]
                self.__vt = self.__wt + self.pos.u

            if "v_r" in input_codes or "w_r" in input_codes:

                if "w_r" in input_codes:

                    i = input_codes.index("w_r")

                else:

                    i = input_codes.index("v_r")

                self.__vr = input_values[i]
                self.__wr = self.__vr
                self.__v = np.sqrt(self.__vt ** 2 + self.__vr ** 2)
                self.__w = np.sqrt(self.__wt ** 2 + self.__wr ** 2)

            elif "w" in input_codes:

                i = input_codes.index("w")
                self.__w = input_values[i]
                self.__wr = np.sqrt(self.__w ** 2 - self.__wt ** 2)
                self.__vr = self.__wr
                self.__v = np.sqrt(self.__vt ** 2 + self.__vr ** 2)

            elif "alpha" in input_codes:

                    i = input_codes.index("alpha")
                    self.__alpha = input_values[i]
                    self.__v = self.__vt / np.sin(self.__alpha)

            elif "beta" in input_codes:

                    i = input_codes.index("beta")
                    self.__beta = input_values[i]
                    self.__wr = self.__wt / np.tan(self.__beta)
                    self.__vr = self.__wr
                    self.__v = np.sqrt(self.__vt ** 2 + self.__vr ** 2)

        elif "v_r" in input_codes or "w_r" in input_codes:

            if "w_r" in input_codes:

                i = input_codes.index("w_r")

            else:

                i = input_codes.index("v_r")

            self.__vr = input_values[i]
            self.__wr = self.__vr

            if "w" in input_codes:

                i = input_codes.index("w")
                self.__w = input_values[i]
                self.__wt = np.sqrt(self.__w ** 2 - self.__wr ** 2)
                self.__vt = self.__wt + self.pos.u
                self.__v = np.sqrt(self.__vt ** 2 + self.__vr ** 2)

            elif "alpha" in input_codes:

                i = input_codes.index("alpha")
                self.__alpha = input_values[i]
                self.__v = self.__vr / np.cos(self.__alpha)
                self.__vt = self.__vr * np.tan(self.__alpha)

            elif "beta" in input_codes:

                i = input_codes.index("beta")
                self.__beta = input_values[i]
                self.__wt = self.__wr * np.tan(self.__beta)
                self.__vt = self.__wt + self.pos.u
                self.__v = np.sqrt(self.__vt ** 2 + self.__vr ** 2)

        elif "w" in input_codes:

            i = input_codes.index("w")
            self.__w = input_values[i]

            if "beta" in input_codes:

                i = input_codes.index("beta")
                self.__beta = input_values[i]
                self.__wt = self.__w * np.sin(self.__beta)
                self.__wr = self.__w * np.cos(self.__beta)
                self.__vt = self.__wt + self.pos.u
                self.__vr = self.__wr
                self.__v = np.sqrt(self.__vt ** 2 + self.__vr ** 2)

    def __evaluate_other(self):

        self.__vr = np.sqrt(self.__v ** 2 - self.__vt ** 2)

        self.__wr = self.__vr
        self.__wt = self.__vt - self.pos.u
        self.__w = np.sqrt(self.__wr ** 2 + self.__wt ** 2)

        self.__beta = np.arcsin(self.__wt / self.w)
        self.__alpha = np.arcsin(self.__vt / self.v)

    @staticmethod
    def __code_is_angle(code):

        return "alpha" == code or "beta" == code

    @staticmethod
    def __code_is_radial(code):

        return "_r" in code

    @staticmethod
    def __code_is_tangential(code):

        return "_t" in code

    def __codes_are_independent(self, code, other_code):

        if self.__code_is_radial(code) and self.__code_is_radial(other_code):
            return False

        if self.__code_is_tangential(code) and self.__code_is_tangential(other_code):
            return False

        if self.__code_is_angle(code) and self.__code_is_angle(other_code):
            return False

        if (code == "alpha" and other_code == "w") or (code == "beta" and other_code == "v"):

            return False

        return self.__check_code_correct(code) and self.__check_code_correct(other_code)
