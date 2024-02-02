from REFPROPConnector import ThermodynamicPoint as TP
from .base_classes.support import BaseTeslaGeometry, BaseTeslaOptions, StatorGeometry, RotorGeometry, StatorOptions, RotorOptions
from .sub_classes.single_phase import Stator, Rotor

""" 

    TODO: 
    
    This python file has been kept to ensure portability with previous versions
    It has to be removed as soon as possible! 
    
"""

class TeslaGeometry(BaseTeslaGeometry):
    def __init__(self):
        super().__init__(stator_geometry=StatorGeometry, rotor_geometry=RotorGeometry)


class TeslaOptions(BaseTeslaOptions):
    def __init__(self):
        super().__init__(stator_options=StatorOptions, rotor_options=RotorOptions)


class TeslaTurbine:

    def __init__(self, fluid, geometry: BaseTeslaGeometry, options: BaseTeslaOptions):

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
