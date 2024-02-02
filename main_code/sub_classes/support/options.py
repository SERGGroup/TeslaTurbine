from main_code.sub_classes.multi_phase.support.void_fraction_handler import void_fraction_options

@void_fraction_options
class RotorOptions:

    profile_rotor = True

    n_rotor = 250
    integr_method = "Std"
    integr_variable = 0.001



class StatorOptions:

    iterate_phi = True
    n_phi_iteration = 50


class TeslaOptions:

    rotor = RotorOptions()
    stator = StatorOptions()
