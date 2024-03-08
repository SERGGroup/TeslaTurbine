from main_code.base_classes.support.options import RotorOptions, StatorOptions, BaseTeslaOptions


class SPRotorOptions(RotorOptions):
    pass

class SPStatorOptions(StatorOptions):
    pass

class SPTeslaOptions(BaseTeslaOptions):

    def __init__(self):

        super().__init__(rotor_options=SPRotorOptions, stator_options=SPStatorOptions)