import matplotlib.pyplot as plt
import numpy as np

data = np.fromregex('fort.10', r"\s*([^\s]+)\s+([^\s]+)\s+([^\s]+)", [('x', float), ('y', float), ('z', float)])
x = data['x']
y = data['y']
z = data['z']

ax = plt.figure(figsize=(20,5)).add_subplot(projection='3d')

ax.scatter(x,y,z)
ax.grid()

plt.show()
