import warnings

DEFAULT_WARNING_MESSAGE = """
   
    !! WARNING from OPTIONS !!
    \"{}\" is not an allowed model.
    Allowed Models are: {}
    The following model has been considered instead: \"{}\"

"""

def append_general_option_modification(

        option_class, option_name,
        acceptable_values, warning_message=DEFAULT_WARNING_MESSAGE,
        default_value=None

):

    if default_value is None:

        default_value = acceptable_values[0]

    storage_name = "__" + option_name

    def models_name_property():

        """Return a property that stores values under a private non-public name."""

        @property
        def prop(self):
            return getattr(self, storage_name)

        @prop.setter
        def prop(self, model):

            model = model.lower()
            if model in acceptable_values:
                setattr(self, storage_name, model)

            else:

                warnings.warn(

                    warning_message.format(

                        model,
                        acceptable_values,
                        getattr(self, storage_name)

                    )

                )

        return prop

    setattr(option_class, storage_name, default_value)
    setattr(option_class, option_name, models_name_property())

    return option_class