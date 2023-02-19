import matplotlib.pyplot as plt
import matplotlib
font = {'weight' : 'bold',
        'size'   : 8}
matplotlib.rc('font', **font)

import numpy as np

data = np.fromregex('fort.2', r"\s*([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)", [('iteration', int), ('kinetic', float), ('potential', float), ('total', float), ('temperature', float)])
kinetic = data['kinetic']
potential = data['potential']
total = data['total']

plt.figure()

dt = 0.004
x = dt * np.arange(1,len(kinetic)+1)

plt.plot(x, kinetic, label='Kinetic')
plt.plot(x, potential, label='Potential')
plt.plot(x, total, label='Total')

plt.xlabel("Time")
plt.grid()
plt.title('Energy')
plt.legend()

plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
            hspace = 0, wspace = 0)
plt.margins(0,0)

plt.savefig(f"energy.pdf", bbox_inches = 'tight', pad_inches = 0, dpi=300)
