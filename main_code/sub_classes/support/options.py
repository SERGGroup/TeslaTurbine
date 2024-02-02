import warnings


class RotorOptions:

    profile_rotor = False
    __TP_EPSILON_MODELS = ["milazzo", "sarti"]

    n_rotor = 250
    integr_method = "Std"
    integr_variable = 0.001
    __tp_epsilon_model = "milazzo"

    @property
    def tp_epsilon_model(self):

        return self.__tp_epsilon_model

    @tp_epsilon_model.setter
    def tp_epsilon_model(self, model: str):

        model = model.lower()
        if model in self.__TP_EPSILON_MODELS:

            self.__tp_epsilon_model = model

        else:

            warnings.warn(

                """
                    !! WARNING from ROTOR OPTIONS !!
                    \"{}\" is not an allowed void fraction model.
                    Allowed Models are: {}
                    The following model has been considered instead: \"{}\"
                
                """.format(

                    model,
                    self.__TP_EPSILON_MODELS,
                    self.__tp_epsilon_model

                )

            )


class StatorOptions:

    iterate_phi = True
    n_phi_iteration = 50


class TeslaOptions:

    rotor = RotorOptions()
    stator = StatorOptions()
