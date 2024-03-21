from REFPROPConnector import ThermodynamicPoint as TP
from .support import BaseTeslaGeometry, BaseTeslaOptions
from .base_stator import BaseStator0D, BaseStator1D
from .base_rotor import BaseRotor


class BaseTeslaTurbine:

    def __init__(

            self, fluid, geometry: BaseTeslaGeometry,
            options: BaseTeslaOptions, stator: type(BaseStator1D),
            rotor: type(BaseRotor)

    ):

        self.fluid = fluid
        self.geometry = geometry
        self.options = options

        self.__init_states(4)

        self.stator = stator(self)
        self.rotor = rotor(self)
        self.P_in = 0.
        self.P_out = 0.
        self.T_in = 0.
    def __init_states(self, n_points):

        base_point = TP([self.fluid], [1], unit_system="MASS BASE SI")

        self.points = list()
        self.static_points = list()

        for i in range(n_points):
            self.points.append(base_point.duplicate())
            self.static_points.append(base_point.duplicate())

    def iterate_pressure(self):
        P_up = self.P_in
        P_down = self.P_out
        n_it = 100

        for i in range (n_it):

            P_1s = (P_up + P_down) / 2

            self.static_points[1].set_variable("P", P_1s)
            self.stator.solve()
            self.rotor.solve()
            P_out_tent = self.static_points[3].get_variable("P")
            error = abs((P_out_tent - self.P_out) / P_out_tent)

            if error < 0.0001:
                break
            elif self.P_out < P_out_tent:
                P_up = P_1s
            else:
                P_down = P_1s

