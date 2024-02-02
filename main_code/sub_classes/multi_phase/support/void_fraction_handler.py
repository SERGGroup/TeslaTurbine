import warnings
import scipy

__ACCEPTABLE_VOID_FRACTION_MODELS = ["milazzo", "sarti"]
__DEFAULT_VOID_FRACTION_MODEL = "milazzo"

def void_fraction_handler(original_class):

    def epsilon(self):

        if self.__epsilon is None:
            self.__epsilon = self.__evaluate_epsilon()

        return self.__epsilon

    def evaluate_epsilon(self):

        x = self.thermo_point.get_variable("x")
        rho_liq = self.liq_phase.get_variable("rho")
        rho_vap = self.vap_phase.get_variable("rho")

        if self.main_rotor.options.tp_epsilon_model == "sarti":

            return 1 / (1 + (1 - x) / x * (rho_vap / rho_liq) ** (2 / 3))

        else:

            m_dot = self.m_dot
            sigma = self.liq_phase.evaluate_RP_code("STN")
            g = scipy.constants.g
            a = (1 + 0.12 * (1 - x)) * (x / rho_vap + (1 - x) / rho_liq)
            b = 1.18 / m_dot * (1 - x) * (g * sigma * (rho_liq - rho_vap) / rho_liq ** 2) ** (1 / 4)

            return x / rho_vap / (a + b)

    def reset_epsilon(self):

        self.__epsilon = None

    setattr(original_class, '__epsilon', None)
    setattr(original_class, 'epsilon', property(epsilon))
    setattr(original_class, 'reset_epsilon', reset_epsilon)
    setattr(original_class, '__evaluate_epsilon', evaluate_epsilon)

    return original_class


def void_fraction_options(option_class):

    storage_name = '__tp_epsilon_model'

    def models_name_property():

        """Return a property that stores values under a private non-public name."""

        @property
        def prop(self):
            return getattr(self, storage_name)

        @prop.setter
        def prop(self, model):

            model = model.lower()
            if model in self.__TP_EPSILON_MODELS:
                setattr(self, storage_name, model)

            else:

                warnings.warn(

                    """
                        !! WARNING from ROTOR OPTIONS !!
                        \"{}\" is not an allowed void fraction model.
                        Allowed Models are: {}
                        The following model has been considered instead: \"{}\"
    
                    """.format(

                        model,
                        __ACCEPTABLE_VOID_FRACTION_MODELS,
                        self.__tp_epsilon_model

                    )

                )

        return prop

    setattr(option_class, storage_name, __DEFAULT_VOID_FRACTION_MODEL)
    setattr(option_class, 'tp_epsilon_model', models_name_property())

    return option_class
