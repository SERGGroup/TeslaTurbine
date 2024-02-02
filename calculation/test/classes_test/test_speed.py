# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.base_classes.support.speed import Position, Speed
import numpy as np

# %%------------   DEFINE TEST PARAM                      -----------------------------------------------------------> #
r = 0.01858
r_new = 0.0175
u = 5.837

v = 28.39
vt = 26.61
vr = 9.906

w = 23.01
wt = 20.77
wr = 9.906

alpha = 69.58
beta = 64.5

theta_0 = 1191
theta_1 = 1200
gamma_0 = 403.7
gamma_1 = 410.8

# %%------------   INIT CLASSES                           -----------------------------------------------------------> #
pos = Position(r, 314.2, theta=theta_0 / 180 * np.pi, gamma=gamma_0 / 180 * np.pi)
spd = Speed(pos)

# %%------------   TEST CLASSES                           -----------------------------------------------------------> #
spd.init_from_codes("v", v, "alpha", alpha / 180 * np.pi)

print("v:     {:6.1f}m/s - {:6.1f}m/s".format(v, spd.v))
print("vr:    {:6.1f}m/s - {:6.1f}m/s".format(vr, spd.vr))
print("vt:    {:6.1f}m/s - {:6.1f}m/s".format(vt, spd.vt))
print("alpha: {:6.1f}°   - {:6.1f}°".format(alpha, spd.alpha / np.pi * 180))
print()
print("w:     {:6.1f}m/s - {:6.1f}m/s".format(w, spd.w))
print("wr:    {:6.1f}m/s - {:6.1f}m/s".format(wr, spd.wr))
print("wt:    {:6.1f}m/s - {:6.1f}m/s".format(wt, spd.wt))
print("beta:  {:6.1f}°   - {:6.1f}°".format(beta, spd.beta / np.pi * 180))

new_pos = spd.get_new_position(r - r_new)
print()
print("t:      {:6.3f} ms".format(new_pos.t * 1000))
print("theta:  {:6.1f}° - {:6.1f}°".format(theta_1, new_pos.theta / np.pi * 180))
print("gamma:  {:6.1f}° - {:6.1f}°".format(gamma_1, new_pos.gamma / np.pi * 180))
