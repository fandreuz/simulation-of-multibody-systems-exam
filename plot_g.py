import matplotlib.pyplot as plt
import matplotlib
font = {'weight' : 'bold',
        'size'   : 20}
matplotlib.rc('font', **font)

import numpy as np

def extract_data(filename):
    return np.fromregex(filename, r"\s*([^\s]+)\s+([^\s]+)", [('rad', float), ('g', float)])

files = [
    # (filename, label)
]

plt.figure()

for filename, label in files:
    data = extract_data(filename)
    r = data['rad']
    g = data['g']

    r = r.reshape((100,-1), order="F")
    g = g.reshape((100,-1), order="F")
    plt.plot(r[:,0], g[:,-1], label=f"T={label}")

plt.legend()
plt.grid()
plt.title("g(r)")

plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
            hspace = 0, wspace = 0)
plt.margins(0,0)

plt.savefig("g.pdf", bbox_inches = 'tight', pad_inches = 0, dpi=300)