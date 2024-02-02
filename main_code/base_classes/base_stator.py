from .support import Position, Speed

class BaseStator:

    def __init__(self, main_turbine):

        self.main_turbine = main_turbine
        self.options = self.main_turbine.options.stator
        self.geometry = self.main_turbine.geometry.stator

        self.input_point = self.main_turbine.points[0]
        self.output_point = self.main_turbine.points[1]

        self.static_input_point = self.main_turbine.static_points[0]
        self.static_output_point = self.main_turbine.static_points[1]
        self.p_out = self.static_output_point.get_variable("P")

        pos = Position(self.geometry.r_int, 0)
        self.speed_out = Speed(pos)

        self.Ma_1 = 0.
        self.m_dot_s = 0.
        self.eta_stat = 0.

    def solve(self):

        pass