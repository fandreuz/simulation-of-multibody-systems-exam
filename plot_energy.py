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

fig = plt.figure()
plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
                    wspace=0.25,
                    hspace=0.25)

dt = 0.004
x = dt * np.arange(1,len(kinetic)+1)

plt.subplot(2,1,1)
plt.plot(x, kinetic)

plt.title('Kinetic energy')

plt.grid()

plt.subplot(2,1,2)
plt.plot(x, potential, label='Potential')
plt.plot(x, total, label='Total')

plt.title('Potential/total energy')

plt.xlabel("Time")
plt.grid()
plt.legend()

plt.savefig(f"energy.pdf", bbox_inches = 'tight', pad_inches = 0, dpi=300)
