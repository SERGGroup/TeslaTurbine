from REFPROPConnector import ThermodynamicPoint as TP
from .support import BaseTeslaGeometry, BaseTeslaOptions
from .base_stator import BaseStator0D, BaseStator1D
from .base_rotor import BaseRotor
import numpy as np


class BaseTeslaTurbine:

    __n_packs = 50
    __m_dot_tot = None

    def __init__(

            self, fluid, geometry: BaseTeslaGeometry,
            options: BaseTeslaOptions, stator: type(BaseStator0D),
            rotor: type(BaseRotor), comp=None

    ):

        if comp is None:

            if type(fluid) == list:
                fluid = fluid[0]

            self.comp = [1]
            self.fluid = [fluid]

        else:

            if type(fluid) == list:

                self.fluid = fluid
                self.comp = comp

            else:

                self.comp = [1]
                self.fluid = [fluid]

        self.geometry = geometry
        self.options = options

        self.__init_states(4)

        self.stator = stator(self)
        self.rotor = rotor(self)

        self.P_in = 0.
        self.P_out = 0.
        self.T_in = 0.

        self.Eta_tesla_ss = 0.
        self.Eta_rotor_ss = 0.
        self.Eta_rotor_ss2 = 0.
        self.work = 0.
        self.power = 0.

    def __init_states(self, n_points):

        base_point = TP(self.fluid, self.comp, unit_system="MASS BASE SI")

        self.points = list()
        self.static_points = list()

        for i in range(n_points):
            self.points.append(base_point.duplicate())
            self.static_points.append(base_point.duplicate())

        self.isentropic_total_outlet = base_point.duplicate()
        self.isentropic_static_outlet = base_point.duplicate()

    def adjusting_flowrate(self, m):

        P_up = self.P_in
        P_down = 200000
        n_it = 100

        for i in range (n_it):

            Ps = (P_up + P_down) / 2

            self.static_points[1].set_variable('P', Ps)

            self.stator.solve()

            m_out = self.stator.m_dot_s
            error = abs((m_out - m) / m)

            if error < 0.00001:
                break
            elif m_out < m:
                P_up = Ps
            else:
                P_down = Ps

        self.rotor.solve()

    def iterate_pressure(self):

        P_up = self.P_in
        P_down = self.P_out
        n_it = 100

        k = 0

        for i in range (n_it):

            P_1s = (P_up + P_down) / 2

            self.solve_with_stator_outlet_pressure(P_1s)

            P_out_tent = self.static_points[3].get_variable("P")
            error = abs((P_out_tent - self.P_out) / P_out_tent)

            if error < 0.0001:
                break
            elif self.P_out < P_out_tent:
                P_up = P_1s
            else:
                P_down = P_1s

            k += 1

            if k == n_it:

                print('Pressure Iteration is not converging')

        # if self.stator.last_stator_speed is not self.stator.out_speed:
        #     print('The stator is CHOKED')

    def iterate_total_pressure(self):

        P_up = self.P_in
        P_down = self.P_out
        n_it = 100

        for i in range (n_it):

            P_1s = (P_up + P_down) / 2

            self.solve_with_stator_outlet_pressure(P_1s)

            P_out_tent = self.points[3].get_variable("P")
            error = abs((P_out_tent - self.P_out) / P_out_tent)

            if error < 0.0001:
                break
            elif self.P_out < P_out_tent:
                P_up = P_1s
            else:
                P_down = P_1s

    def solve_with_stator_outlet_pressure(self, P_1s):

        self.static_points[1].set_variable("P", P_1s)
        self.stator.solve()
        self.rotor.solve()

    def fix_rotor_inlet_condition(self, P_1s, T_1s, alpha, v):

        # self.static_points[1].set_variable("P", P_1s)
        # self.stator.solve()

        self.stator.speed_out.init_from_codes("alpha", alpha, "v", v)
        self.static_points[1].set_variable("P", P_1s)
        self.static_points[1].set_variable("T", T_1s)

        self.points[1].set_variable("h", self.static_points[1].get_variable("h") + v ** 2 / 2)
        self.points[1].set_variable("s", self.static_points[1].get_variable("s"))

        self.rotor.solve()

    def fix_rotor_inlet_condition2(self, P_1s, T_1s, u, v):

        self.stator.speed_out.init_from_codes("v_r", u, "v", v)
        self.static_points[1].set_variable("P", P_1s)
        self.static_points[1].set_variable("T", T_1s)

        self.points[1].set_variable("h", self.static_points[1].get_variable("h") + v ** 2 / 2)
        self.points[1].set_variable("s", self.static_points[1].get_variable("s"))

        self.rotor.solve()

    def fix_rotor_inlet_condition_quality(self, P_1s, x_1s, alpha, v):

        # self.static_points[1].set_variable("P", P_1s)
        # self.stator.solve()

        self.stator.speed_out.init_from_codes("alpha", alpha, "v", v)
        self.static_points[1].set_variable("P", P_1s)
        self.static_points[1].set_variable("x", x_1s)

        self.points[1].set_variable("h", self.static_points[1].get_variable("h") + v ** 2 / 2)
        self.points[1].set_variable("s", self.static_points[1].get_variable("s"))
        # self.points[1].set_variable("x", 0)

        self.rotor.solve()

    def evaluate_performances(self):

        self.static_points[0].copy_state_to(self.isentropic_total_outlet)
        self.static_points[0].copy_state_to(self.isentropic_static_outlet)

        self.isentropic_total_outlet.set_variable("s", self.points[2].get_variable("s"))
        self.isentropic_total_outlet.set_variable("p", self.points[3].get_variable("p"))

        self.isentropic_static_outlet.set_variable("s", self.static_points[2].get_variable("s"))
        self.isentropic_static_outlet.set_variable("p", self.static_points[3].get_variable("p"))

        self.Eta_tesla_ts = (self.points[0].get_variable("h")-self.static_points[3].get_variable("h")) / (
                self.points[0].get_variable("h") - self.isentropic_static_outlet.get_variable("h"))

        self.Eta_rotor_tt = (self.points[0].get_variable("h") - self.points[3].get_variable("h")) / (
                self.points[0].get_variable("h") - self.isentropic_total_outlet.get_variable("h"))

        self.Eta_tesla_ss = (self.static_points[0].get_variable("h")-self.static_points[3].get_variable("h")) / (
                self.static_points[0].get_variable("h") - self.isentropic_static_outlet.get_variable("h"))

        self.work = self.rotor.first_speed.vt * self.rotor.first_speed.u - self.rotor.rotor_points[-1].speed.vt * self.rotor.rotor_points[-1].speed.u
        self.work2 = self.points[2].get_variable("h") - self.points[3].get_variable("h")
        self.power = self.work * self.rotor.m_dot_ch

        self.Eta_tesla_tt2 = self.work / (self.points[0].get_variable("h") - self.isentropic_total_outlet.get_variable("h"))

    def evaluate_bearing_performances(self, M_lost):

        self.static_points[0].copy_state_to(self.isentropic_outlet)

        self.isentropic_outlet.set_variable("s", self.static_points[0].get_variable("s"))
        self.isentropic_outlet.set_variable("p", self.static_points[3].get_variable("p"))

        self.Eta_tesla_ts = (self.points[0].get_variable("h")-self.points[3].get_variable("h")) / (
                self.static_points[0].get_variable("h") - self.isentropic_outlet.get_variable("h"))

        self.Eta_rotor_tt = (self.static_points[0].get_variable("h") - self.static_points[3].get_variable("h")) / (
                self.static_points[0].get_variable("h") - self.isentropic_outlet.get_variable("h"))

        self.Eta_tesla_ss = (self.static_points[0].get_variable("h")-self.static_points[3].get_variable("h")) / (
                self.static_points[0].get_variable("h") - self.isentropic_outlet.get_variable("h"))

        self.work = self.rotor.first_speed.vt * self.rotor.first_speed.u - self.rotor.rotor_points[-1].speed.vt * self.rotor.rotor_points[-1].speed.u
        self.power = self.work * self.rotor.m_dot_ch

        self.Eta_tesla_ss = self.work / (self.static_points[0].get_variable("h") - self.isentropic_outlet.get_variable("h"))

        self.Eta_bearings = (self.work - M_lost * self.rotor.omega / (self.rotor.m_dot_ch * self.geometry.rotor.n_discs * self.n_packs)) / (self.static_points[0].get_variable("h") - self.isentropic_outlet.get_variable("h"))

    def evaluate_rotor_performances(self):

        self.isentropic_outlet = self.static_points[2].duplicate()

        self.isentropic_outlet.set_variable("s", self.points[2].get_variable("s"))
        self.isentropic_outlet.set_variable("p", self.points[3].get_variable("p"))

        # self.Eta_rotor_ss = (self.static_points[2].get_variable("h") - self.static_points[3].get_variable("h")) / (
        #         self.static_points[2].get_variable("h") - self.isentropic_outlet.get_variable("h"))
        #
        # self.work = self.rotor.first_speed.vt * self.rotor.first_speed.u - self.rotor.rotor_points[-1].speed.vt * \
        #             self.rotor.rotor_points[-1].speed.u
        # self.power = self.work * self.rotor.m_dot_ch
        self.work = self.rotor.first_speed.vt * self.rotor.first_speed.u - self.rotor.rotor_points[-1].speed.vt * self.rotor.rotor_points[-1].speed.u
        self.work2 = self.points[2].get_variable("h")-self.points[3].get_variable("h")
        self.power = self.work * self.rotor.m_dot_ch

        self.Eta_rotor_tt = self.work / (self.points[2].get_variable("h") - self.isentropic_outlet.get_variable("h"))

    @property
    def m_dot_tot(self):

        return self.rotor.m_dot_ch * self.geometry.rotor.n_discs * self.__n_packs

    @m_dot_tot.setter
    def m_dot_tot(self, m_dot_tot):

        self.__m_dot_tot = m_dot_tot
        self.__n_packs = None

    @property
    def n_packs(self):

        if self.__m_dot_tot is not None:

            try:

                # round up to the highest integer
                # return np.ceil(self.__m_dot_tot / (self.rotor.m_dot_ch * self.n_packs))
                return np.ceil(self.__m_dot_tot / self.rotor.m_dot_ch)

            except:

                return None

        else:

            return self.__n_packs

    @n_packs.setter
    def n_packs(self, n_packs):

        self.__m_dot_tot = None
        self.__n_packs = n_packs

    @property
    def volume(self):

        return np.pi * self.geometry.stator.d_0 ** 2 / 4 * (self.geometry.H_s * self.n_packs)