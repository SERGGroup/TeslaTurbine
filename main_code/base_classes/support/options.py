class RotorOptions:

    profile_rotor = True
    n_rotor = 250
    integr_method = "Std"
    integr_variable = 0.001


class StatorOptions:

    iterate_phi = True
    n_phi_iteration = 50

    profile_stator = True
    n_stator = 101
    integr_method = "Std"
    roughness = 5e-7


class BaseTeslaOptions:

    def __init__(

            self, rotor_options: type(RotorOptions),
            stator_options: type(StatorOptions)

    ):

        self.rotor = rotor_options()
        self.stator = stator_options()
