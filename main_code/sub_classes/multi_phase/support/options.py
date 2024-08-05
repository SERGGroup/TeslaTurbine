from main_code.base_classes.support.options import RotorOptions, StatorOptions, BaseTeslaOptions
from main_code.sub_classes.multi_phase.handlers.void_fraction_handler import void_fraction_options

@void_fraction_options
class TPRotorOptions(RotorOptions):
    pass

class TPStatorOptions(StatorOptions):
    pass

class TPTeslaOptions(BaseTeslaOptions):

    def __init__(self):

        super().__init__(rotor_options=TPRotorOptions, stator_options=TPStatorOptions)
