import numpy as np
from REFPROPConnector import ThermodynamicPoint as TP

class TeslaGeometry:

    throat_width = 0.0002
    d_rotor = 2000

class TeslaOptions:

    n_rotor = 250
    profile_rotor = False

class TeslaTurbine:

    def __init__(self, fluid, geometry: TeslaGeometry, options: TeslaOptions):

        self.fluid = fluid
        self.geometry = geometry
        self.options = options

        self.__init_states(4)
        self.__init_profilers()

    def __init_states(self, n_points):

        base_point = TP([self.fluid], [1])

        self.points = list()
        self.static_points = list()

        for i in range(n_points):
            self.points.append(base_point.duplicate())
            self.static_points.append(base_point.duplicate())

    def __init_profilers(self):

        self.rotor_profiler = np.array([])

    def iterate_pressure(self):

        pass

    def stator(self):

        self.points[0].set_variable("T", 10)
        self.points[0].set_variable("P", 1)

        print(self.geometry.throat_width)
        print(self.points[0].get_unit("H"))
        print(self.points[0].evaluate_RP_code("D2DDPT"))



if __name__ == "__main__":

    curr_geometry = TeslaGeometry()
    curr_options = TeslaOptions()
    tt = TeslaTurbine("Carbon Dioxide", curr_geometry, curr_options)

    for throat in range(5):

        curr_geometry.throat_width = throat
        tt.stator()
