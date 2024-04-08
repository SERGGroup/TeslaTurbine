from .base_rotor import BaseRotor, BaseRotorStep
from .base_stator import BaseStator0D, BaseStatorStep, BaseStator1D
from .base_turbine import BaseTeslaTurbine
from .support import (

    BaseTeslaGeometry, BaseTeslaOptions,
    StatorOptions, StatorGeometry,
    RotorGeometry, RotorOptions,
    Position, Speed

)

import os.path
CODE_FOLDER = os.path.dirname(os.path.dirname(__file__))
CALCULATION_FOLDER = os.path.join(os.path.dirname(CODE_FOLDER), 'calculation')
