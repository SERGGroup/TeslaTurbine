class RotorOptions:

    n_rotor = 250
    profile_rotor = False

    integr_method = "Std"
    integr_variable = 2


class StatorOptions:

    iterate_phi = True
    n_phi_iteration = 50


class TeslaOptions:

    rotor = RotorOptions()
    stator = StatorOptions()
