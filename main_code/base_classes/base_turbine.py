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

    def __init_states(self, n_points):

        base_point = TP([self.fluid], [1], unit_system="MASS BASE SI")

        self.points = list()
        self.static_points = list()

        for i in range(n_points):
            self.points.append(base_point.duplicate())
            self.static_points.append(base_point.duplicate())

    def iterate_pressure(self):

        # TODO
        pass
