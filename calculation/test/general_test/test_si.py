# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
from REFPROPConnector import ThermodynamicPoint as TP

# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
base_point = TP(["CarbonDioxide"], [1], unit_system="MASS BASE SI")
base_point.list_properties()

# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
new_point = TP(["CarbonDioxide"], [1])
new_point.list_properties()

# %%------------   IMPORT CLASSES                         -----------------------------------------------------------> #
tk = 300
base_point.set_variable("T", tk)
base_point.set_variable("P", 101000)

base_point.copy_state_to(new_point)

print(new_point.get_variable("P"))
print(new_point.get_variable("T"))

