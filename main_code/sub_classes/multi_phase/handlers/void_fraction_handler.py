from main_code.utils.general_option_modifier import append_general_option_modification
import numpy as np

__ACCEPTABLE_VOID_FRACTION_MODELS = ["chisholm", "sarti"]
__DEFAULT_VOID_FRACTION_MODEL = "chisholm"
__OPTION_NAME = 'tp_epsilon_model'

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

    Inside the class there should be:

        - a "self.m_dot" variable: A float containing the flow rate of the fluid
        - a "self.thermo_point" variable: A REFPROP_connector "ThermodynamicPoint" containing the thermodynamic state.
        - a "self.liq_phase" ad a "self.vap_phase" variable: two REFPROP_connector "ThermodynamicPoint" containing the
        liquid and vapor properties of the fluid.

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
        rho_ratio = rho_vap / rho_liq
        if getattr(self.options, __OPTION_NAME) == "sarti":

            """
                
                From:
                
                    G.Sarti, "Development of an experimental data reduction method and 
                    a two-phase model of a Tesla turbine for organic fluid", Master Thesis, 2020
                    
                Could be found in:
                
                    documentation/literature/Multi-Phase/Thesis/Tesi Sarti/Final Material/Tesi_Sarti_FullText.pdf
            
            """
            return 1 / (1 + (1 - x) / x * rho_ratio ** (2 / 3))

        else:

            """

                From:

                    D. Chisholm, "RESEARCH NOTE: VOID FRACTION DURING TWO-PHASE FLOW", 1973

                Could be found in:

                    documentation/literature/Multi-Phase/Void Fraction Modelling/chisholm1973 -other.pdf
                    
            """

            return 1 / (1 + (1 - x) / x * rho_ratio * np.sqrt(1 - x * (1 - 1 / rho_ratio)))

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

    warning_message = """
    
        !! WARNING from OPTIONS !!
        \"{}\" is not an allowed void fraction model.
        Allowed Models are: {}
        The following model has been considered instead: \"{}\"

    """
    option_class = append_general_option_modification(

        option_class=option_class, option_name=__OPTION_NAME,
        acceptable_values=__ACCEPTABLE_VOID_FRACTION_MODELS,
        warning_message=warning_message, default_value=__DEFAULT_VOID_FRACTION_MODEL

    )

    return option_class
