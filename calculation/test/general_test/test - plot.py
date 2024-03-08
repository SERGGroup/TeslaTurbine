# %%------------   IMPORT CLASSES

import numpy as np
import matplotlib.pyplot as plt

# %%------------   PLOT PRESSURE LOSSES

BarWidth = 0.2
fig = plt.subplots(figsize=(10, 6))

X1 = [9.099, 17.564, 21.974, 30.089, 19.991, 19.371]
X2 = [34.046, 66.058, 82.851, 119, 79.065, 76.612]

br1 = np.arange(len(X1))
br2 = [x + BarWidth + 0.02 for x in br1]

plt.bar(br1, X1, color ='darkblue', width = BarWidth, edgecolor = 'darkblue', label = 'Supply')
plt.bar(br2, X2, color = 'darkorange', width = BarWidth, edgecolor = 'darkorange', label = 'Return')

plt.xlabel('Section', fontweight = 'bold', fontsize = 15)
plt.ylabel('Pressure Losses [kPa]', fontweight = 'bold', fontsize = 15)
plt.xticks([r + BarWidth/2 for r in range(len(X1))], ["1","2", "3", "4", "5", "6"], fontsize = '12')
plt.yticks(fontsize = '12')

plt.legend(loc = 'upper right', fontsize = '15', edgecolor = 'black')
plt.show()

# %%------------   PLOT COST DISTRIBUTION (WATER)

costs_h2o = [8980000, 1360000, 6193, 1376000]
labels_h2o = ["Well", "Pump", "Pipeline", "Heat Exchangers"]

fig1, ax1 = plt.subplots()
ax1.pie(costs_h2o, labels = labels_h2o, colors = ['saddlebrown', 'peru', 'darkgrey', 'firebrick'], autopct='%1.1f%%')

plt.show()

# %%------------   PLOT COST DISTRIBUTION (WATER)

costs_co2 = [4500000, 1714000, 11632, 5928000]
labels_co2 = ["Well", "Heat Pump", "Pipeline", "Heat Exchangers"]

fig2, ax2 = plt.subplots()
ax2.pie(costs_co2, labels = labels_co2, colors = ['saddlebrown', 'peru', 'darkgrey', 'firebrick'], autopct='%1.1f%%')

plt.show()

# %%------------   PLOT COST DISTRIBUTION (COMPARISON)

BarWidth1 = 0.2
fig3 = plt.subplots(figsize=(10, 6))

X3 = [1360, 8980, 1376,  11720]
X4 = [1714, 4500, 5928, 12150]

br3 = np.arange(len(X3))
br4 = [x + BarWidth1 + 0.02 for x in br3]

plt.bar(br3, X3, color ='darkblue', width = BarWidth1, edgecolor = 'darkblue', label = 'Water')
plt.bar(br4, X4, color = 'darkorange', width = BarWidth1, edgecolor = 'darkorange', label = 'CO2')

# plt.xlabel('Component', fontweight = 'bold', fontsize = 15)
plt.ylabel('Cost [k€]', fontweight = 'bold', fontsize = 15)
plt.xticks([r + BarWidth1/2 for r in range(len(X3))], ["Pumps/Compressors", "Well", "Heat Exchangers", "Total"], fontsize = '15', fontweight = 'bold')
plt.yticks(fontsize = '12')

plt.legend(loc = 'upper left', fontsize = '15', edgecolor = 'black')
plt.show()