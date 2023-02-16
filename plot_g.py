import matplotlib.pyplot as plt
import numpy as np

data = np.fromregex('fort.9', r"\s*([^\s]+)\s+([^\s]+)", [('rad', float), ('g', float)])
r = data['rad']
g = data['g']
plt.figure(figsize=(20,5))

plt.plot(r,g)
plt.grid()

plt.title('g(r)')

plt.show()
