# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from main_code.support.speed_calculator import Position, Speed
import numpy as np

# %%------------   DEFINE TEST PARAM                      -----------------------------------------------------------> #
r = 0.01858
u = 5.837

v = 28.39
vt = 26.61
vr = 9.906

w = 23.01
wt = 20.77
wr = 9.906

alpha = 69.58
beta = 64.5

# %%------------   INIT CLASSES                           -----------------------------------------------------------> #
pos = Position(r, 0, 0, 314.2)
spd = Speed(pos)

# %%------------   TEST CLASSES                           -----------------------------------------------------------> #
spd.init_from_codes("v", v, "w", w)

print("v:     {:.3f} - {:.3f}".format(v, spd.v))
print("vr:    {:.3f} - {:.3f}".format(vr, spd.vr))
print("vt:    {:.3f} - {:.3f}".format(vt, spd.vt))
print("alpha: {:.3f} - {:.3f}".format(alpha, spd.alpha / np.pi * 180))
print()
print("w:     {:.3f} - {:.3f}".format(w, spd.w))
print("wr:    {:.3f} - {:.3f}".format(wr, spd.wr))
print("wt:    {:.3f} - {:.3f}".format(wt, spd.wt))
print("beta:  {:.3f} - {:.3f}".format(beta, spd.beta / np.pi * 180))
