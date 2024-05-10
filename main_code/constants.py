import os

CURRENT_FOLDER = os.path.dirname(os.path.abspath(__file__))
MAIN_FOLDER = os.path.abspath(os.path.join(CURRENT_FOLDER, os.pardir))
CALCULATION_FOLDER = os.path.join(MAIN_FOLDER, "calculation")