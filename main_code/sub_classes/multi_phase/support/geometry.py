from main_code.base_classes.support.geometry import BaseTeslaGeometry, RotorGeometry, StatorGeometry


class TPRotorGeometry(RotorGeometry):
    pass


class TPStatorGeometry(StatorGeometry):
    pass


class TPTeslaGeometry(BaseTeslaGeometry):

    def __init__(self):

        super().__init__(rotor_geometry=TPRotorGeometry, stator_geometry=TPStatorGeometry)
