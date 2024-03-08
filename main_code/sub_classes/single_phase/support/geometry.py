from main_code.base_classes.support.geometry import BaseTeslaGeometry, RotorGeometry, StatorGeometry


class SPRotorGeometry(RotorGeometry):
    pass


class SPStatorGeometry(StatorGeometry):
    pass


class SPTeslaGeometry(BaseTeslaGeometry):

    def __init__(self):

        super().__init__(rotor_geometry=SPRotorGeometry, stator_geometry=SPStatorGeometry)