from sub_classes.support import TeslaGeometry, TeslaOptions
from REFPROPConnector import ThermodynamicPoint as TP
from sub_classes.multi_phase import Stator, Rotor
import numpy as np


class TeslaTurbine:

    def __init__(self, fluid, geometry: TeslaGeometry, options: TeslaOptions):

        self.fluid = fluid
        self.geometry = geometry
        self.options = options

        self.stator = Stator(self)
        self.rotor = Rotor(self)

        self.__init_states(4)

    def __init_states(self, n_points):

        base_point = TP([self.fluid], [1])

        self.points = list()
        self.static_points = list()

        for i in range(n_points):
            self.points.append(base_point.duplicate())
            self.static_points.append(base_point.duplicate())

    def iterate_pressure(self):

        pass


if __name__ == "__main__":

    curr_geometry = TeslaGeometry()
    curr_options = TeslaOptions()
    tt = TeslaTurbine("Carbon Dioxide", curr_geometry, curr_options)

    for throat in range(5):

        curr_geometry.stator.throat_width = throat
        tt.stator()

    print(tt.geometry.stator.r_int)
