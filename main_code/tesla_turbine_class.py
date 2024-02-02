from REFPROPConnector import ThermodynamicPoint as TP
from .sub_classes.support import (TeslaGeometry, TeslaOptions)
from .sub_classes.multi_phase import Stator, Rotor
import numpy as np


class TeslaTurbine:

    def __init__(self, fluid, geometry: TeslaGeometry, options: TeslaOptions):

        self.fluid = fluid
        self.geometry = geometry
        self.options = options

        self.__init_states(4)

        self.stator = Stator(self)
        self.rotor = Rotor(self)

    def __init_states(self, n_points):

        base_point = TP([self.fluid], [1], unit_system="MASS BASE SI")

        self.points = list()
        self.static_points = list()

        for i in range(n_points):
            self.points.append(base_point.duplicate())
            self.static_points.append(base_point.duplicate())

    def iterate_pressure(self):

        pass
