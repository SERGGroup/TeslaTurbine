from main_code.base_classes.base_rotor import BaseRotor, BaseRotorStep, Speed

class SPRotorStep(BaseRotorStep):

    def __init__(self, main_rotor, speed: Speed):

        super().__init__(main_rotor, speed)

    def get_variations(self, dr):

        return 0., 0., 0., 0.

    def solve(self):

        pass

class SPRotor(BaseRotor):

    def __init__(self, main_turbine):

        super().__init__(main_turbine, SPRotorStep)