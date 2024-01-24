from REFPROPConnector import ThermodynamicPoint as TP

class TeslaGeometry:

    throath_width = 0.0002


class TeslaTurbine:

    def __init__(self, fluid, geometry: TeslaGeometry):

        self.fluid = fluid
        self.geometry = geometry
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

    def stator(self):

        self.points[0].set_variable("T", 10)
        self.points[0].set_variable("P", 1)

        print(self.geometry.troath_width)
        print(self.points[0].get_unit("H"))
        print(self.points[0].evaluate_RP_code("D2DDPT"))


if __name__ == "__main__":
    geometry = TeslaGeometry()
    tt = TeslaTurbine("Carbon Dioxide", geometry)
    for throat in range(5):
        geometry.throath_width = throat
        tt.stator()

