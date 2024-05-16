from REFPROPConnector import ThermodynamicPoint as TP
from .support import BaseTeslaGeometry, BaseTeslaOptions
from .base_diffuser import BaseDiffuser, NoDiffuser
from .base_stator import BaseStator0D
from .base_rotor import BaseRotor
import numpy as np


class BaseTeslaTurbine:

    __n_packs = 50
    __m_dot_tot = None

    def __init__(

            self, fluid, geometry: BaseTeslaGeometry,
            options: BaseTeslaOptions, stator: type(BaseStator0D),
            rotor: type(BaseRotor), diffuser: type(BaseDiffuser) = NoDiffuser

    ):

        self.fluid = fluid
        self.geometry = geometry
        self.options = options

        self.__init_states(5)

        self.stator = stator(self)
        self.rotor = rotor(self)
        self.diffuser = diffuser(self)

        self.P_in = 0.
        self.P_out = 0.
        self.T_in = 0.

        self.eta_tt = 0.
        self.eta_ss = 0.

        self.work = 0.
        self.power = 0.

    def __init_states(self, n_points):

        base_point = TP([self.fluid], [1], unit_system="MASS BASE SI")

        self.points = list()
        self.static_points = list()

        for i in range(n_points):
            self.points.append(base_point.duplicate())
            self.static_points.append(base_point.duplicate())

        self.iso_entropic_outlet = base_point.duplicate()

    def iterate_pressure(self):

        P_up = self.P_in
        P_down = self.P_out
        n_it = 100

        for i in range (n_it):

            P_1s = (P_up + P_down) / 2

            try:

                self.solve_with_stator_outlet_pressure(P_1s, check_P_out=True)

                P_out_tent = self.static_points[-1].get_variable("P")
                error = abs((P_out_tent - self.P_out) / P_out_tent)

                if error < 0.0001:
                    break
                elif self.P_out < P_out_tent:
                    P_up = P_1s
                else:
                    P_down = P_1s

            except:

                P_down = P_1s

    def solve_with_stator_outlet_pressure(self, P_1s, check_P_out=False):

        self.static_points[1].set_variable("P", P_1s)
        self.stator.solve()
        self.rotor.solve()

        if (not check_P_out) or (self.rotor.output_point.get_variable("P") > self.P_out):

            self.diffuser.solve()

        else:

            self.rotor.output_point.copy_state_to(self.points[-1])
            self.rotor.static_output_point.copy_state_to(self.static_points[-1])

    def fix_rotor_inlet_condition(self, P_1s, T_1s, alpha, v):

        self.static_points[1].set_variable("P", P_1s)
        self.stator.solve()

        self.stator.speed_out.init_from_codes("alpha", alpha, "v", v)
        self.static_points[1].set_variable("P", P_1s)
        self.static_points[1].set_variable("T", T_1s)
        self.rotor.solve()

    def evaluate_performances(self):

        self.static_points[0].copy_state_to(self.iso_entropic_outlet)

        self.iso_entropic_outlet.set_variable("s", self.static_points[0].get_variable("s"))
        self.iso_entropic_outlet.set_variable("p", self.static_points[3].get_variable("p"))

        self.eta_ss = (self.points[0].get_variable("h") - self.points[3].get_variable("h")) / (
                self.static_points[0].get_variable("h") - self.iso_entropic_outlet.get_variable("h"))

        self.work = self.rotor.first_speed.vt * self.rotor.first_speed.u - self.rotor.rotor_points[-1].speed.vt * self.rotor.rotor_points[-1].speed.u
        self.power = self.work * self.rotor.m_dot_ch

        self.eta_tt = self.work / (self.static_points[0].get_variable("h") - self.iso_entropic_outlet.get_variable("h"))

    @property
    def m_dot_tot(self):

        return self.rotor.m_dot_ch * self.geometry.rotor.n_channels * self.__n_packs

    @m_dot_tot.setter
    def m_dot_tot(self, m_dot_tot):

        self.__m_dot_tot = m_dot_tot
        self.__n_packs = None

    @property
    def n_packs(self):

        if self.__m_dot_tot is not None:

            try:

                # round up to the highest integer
                return np.ceil(self.__m_dot_tot / (self.rotor.m_dot_ch * self.geometry.rotor.n_channels))

            except:

                return None

        else:

            return self.__n_packs

    @n_packs.setter
    def n_packs(self, n_packs):

        self.__m_dot_tot = None
        self.__n_packs = n_packs
