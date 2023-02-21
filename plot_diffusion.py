import matplotlib.pyplot as plt
import matplotlib
import numpy as np

font = {'weight' : 'bold',
        'size'   : 20}
matplotlib.rc('font', **font)

def extract_data(filename):
    return np.fromregex(filename, r"\s*([^\s]+)\s+([^\s]+)", [('i', int), ('r2', float)])

def compute_D(r2):
    return r2[-1] / (6 * dt * len(r2))

files = [
    # (filename, x)
]

plt.figure()

dt = 0.004
skip = 500

Ds = []
for filename, _ in files:
    data = extract_data(filename)
    Ds.append(compute_D(data['r2']))

xs = list(zip(*files))[1]
plt.plot(xs, Ds)

plt.xlabel("T0")

plt.grid()
plt.title("$D$")

plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
            hspace = 0, wspace = 0)
plt.margins(0,0)
plt.savefig("diffusion.pdf", bbox_inches = 'tight', pad_inches = 0, dpi=300)