import warnings
import scipy.constants

__ACCEPTABLE_VOID_FRACTION_MODELS = ["milazzo", "sarti"]
__DEFAULT_VOID_FRACTION_MODEL = "milazzo"

"""
    
    This module implement the void fraction calculation inside a class (usually the rotor or stator step) 
    by creating an "epsilon" property for the desired class and by defining the possible ways of performing 
    the evaluation.
    
    The modules also automatically updates the option class by appending a property that allow to chose between 
    the acceptable models. 
    
"""


def void_fraction_handler(original_class):

    """
    This module implement the void fraction calculation inside a class (usually the rotor or stator step).
    This is done by adding the following elements to the class:

        - an "__epsilon" variable: stores the calculated void fraction (by default is None).
        - an "evaluate_epsilon" method: that contains the epsilon evaluation procedure.
        - an "epsilon" property: can be used to retrieve the calculated void fraction. If "__epsilon" is None then
        evaluate it by calling the "evaluate_epsilon()" method.
        - a "reset_epsilon" method: reset "__epsilon" to None so that, once the "epsilon" property is called again,
        it will invoke again the "evaluate_epsilon" method.

    these properties will be added to the class upon class generation using the decorator in the following way:

        @void_fraction_handler

        class <ClassName>:

            ...

    """

    def epsilon(self):

        if self.__epsilon is None:
            self.__epsilon = self.__evaluate_epsilon()

        return self.__epsilon

    def evaluate_epsilon(self):

        x = self.thermo_point.get_variable("x")
        rho_liq = self.liq_phase.get_variable("rho")
        rho_vap = self.vap_phase.get_variable("rho")

        if self.options.tp_epsilon_model == "sarti":

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

    """
    This module add the options related to the void fraction modelling to the option class of the Rotor or the Stator
    This is done by adding the following elements to the class:

        - an "__tp_epsilon_model" variable: stores the name of the model to be used in the calculation.
        - a "tp_epsilon_model" property: to set and get the value of "__tp_epsilon_model". The value is set
        only if the provided model name is part of the "__ACCEPTABLE_VOID_FRACTION_MODELS" list.

    these properties will be added to the class upon class generation using the decorator in the following way:

        @void_fraction_options

        class <OptionClassName>:

            ...

    """

    storage_name = '__tp_epsilon_model'

    def models_name_property():

        """Return a property that stores values under a private non-public name."""

        @property
        def prop(self):
            return getattr(self, storage_name)

        @prop.setter
        def prop(self, model):

            model = model.lower()
            if model in __ACCEPTABLE_VOID_FRACTION_MODELS:
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
