class RotorOptions:

    profile_rotor = False

    n_rotor = 250
    integr_method = "Std"
    integr_variable = 0.001


class StatorOptions:

    iterate_phi = True
    n_phi_iteration = 50


class TeslaOptions:

    rotor = RotorOptions()
    stator = StatorOptions()
