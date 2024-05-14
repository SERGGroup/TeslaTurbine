from main_code.utils.general_option_modifier import append_general_option_modification
import numpy as np

__ACCEPTABLE_MODELS = ["chisholm", "sarti"]
__DEFAULT_MODEL = "chisholm"
__OPTION_NAME = 'flow_losses_model'

"""

    This module implement the pipe losses calculation inside a class (usually the rotor or stator step) 
    by creating a "dpdl" method for the desired class and by defining the possible ways of performing 
    the evaluation.

    The modules also automatically updates the option class by appending a property that allow to chose between 
    the acceptable models. 

"""


def flow_losses_handler(original_class):
    """

    TODO

    these properties will be added to the class upon class generation using the decorator in the following way:

        @void_fraction_handler

        class <ClassName>:

            ...

    """
    def dpdl(self):

        x = self.total_point.get_variable("x")
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

    setattr(original_class, 'dpdl', dpdl)
    return original_class


def flow_losses_options(option_class):
    """
    This module add the options related to the pipe flow losses modelling to the option class of the Rotor or the Stator
    This is done by adding the following elements to the class:

        - an "__flow_losses_model" variable: stores the name of the model to be used in the calculation.
        - a "flow_losses_model" property: to set and get the value of "__flow_losses_model". The value is set
        only if the provided model name is part of the "__ACCEPTABLE_VOID_FRACTION_MODELS" list.

    these properties will be added to the class upon class generation using the decorator in the following way:

        @void_fraction_options

        class <OptionClassName>:

            ...

    """

    warning_message = """

        !! WARNING from OPTIONS !!
        \"{}\" is not an allowed pipe flow losses model.
        Allowed Models are: {}
        The following model has been considered instead: \"{}\"

    """
    option_class = append_general_option_modification(

        option_class=option_class, option_name=__OPTION_NAME,
        acceptable_values=__ACCEPTABLE_MODELS, warning_message=warning_message,
        default_value=__DEFAULT_MODEL

    )

    return option_class

