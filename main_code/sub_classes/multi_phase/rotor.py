from main_code.base_classes.base_rotor import BaseRotor, BaseRotorStep, Speed
from .handlers.flow_pressure_losses_handler import flow_losses_handler
from .handlers.void_fraction_handler import void_fraction_handler

@flow_losses_handler
@void_fraction_handler
class TPRotorStep(BaseRotorStep):

    def __init__(self, main_rotor, speed: Speed):

        super().__init__(main_rotor, speed)

        self.liq_speed = Speed(speed.pos)
        self.vap_speed = Speed(speed.pos)

        self.liq_phase = self.main_rotor.input_point.duplicate()
        self.vap_phase = self.main_rotor.input_point.duplicate()

    def get_variations(self, dr):

        self.evaluate_condition()
        return 0., 0., 0., 0.

    def evaluate_condition(self):

        self.reset_epsilon()

        self.liq_phase.set_variable("P", self.thermo_point.get_variable("P"))
        self.vap_phase.set_variable("P", self.thermo_point.get_variable("P"))

        self.liq_phase.set_variable("x", 0)
        self.vap_phase.set_variable("x", 1)

        print(self.epsilon)


class TPRotor(BaseRotor):

    def __init__(self, main_turbine):

        super().__init__(main_turbine, TPRotorStep)
