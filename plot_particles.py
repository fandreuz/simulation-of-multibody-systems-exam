import matplotlib.pyplot as plt
import matplotlib
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 8}
matplotlib.rc('font', **font)

import numpy as np

data = np.fromregex('fort.10', r"\s*([^\s]+)\s+([^\s]+)\s+([^\s]+)", [('x', float), ('y', float), ('z', float)])
x = data['x']
y = data['y']
z = data['z']

ax = plt.figure().add_subplot(projection='3d')

ax.scatter(x,y,z, edgecolors="black")
ax.grid()

plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
            hspace = 0, wspace = 0)
plt.margins(0,0,0)
plt.savefig("particles_initial.pdf", bbox_inches = 'tight', pad_inches = 0, dpi=300)
