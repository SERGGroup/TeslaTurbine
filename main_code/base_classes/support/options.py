class RotorOptions:

    profile_rotor = True
    n_rotor = 250
    integr_method = "Std"
    integr_variable = 0.001
    sp_check = False


class StatorOptions:

    iterate_phi = True
    n_phi_iteration = 150
    metastability_check = False

    profile_stator = True
    n_stator = 101
    integr_method = "Std"
    roughness = 5e-7

class DiffuserOptions:

    n_max_iter = 20
    alpha = 0.9
    toll = 1e-5

class BaseTeslaOptions:

    def __init__(

            self, rotor_options: type(RotorOptions),
            stator_options: type(StatorOptions),
            diffuser_options: type(DiffuserOptions) = DiffuserOptions,

    ):

        self.rotor = rotor_options()
        self.stator = stator_options()
        self.diffuser = diffuser_options()
