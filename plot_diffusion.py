import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from sklearn.linear_model import LinearRegression

font = {'weight' : 'bold',
        'size'   : 20}
matplotlib.rc('font', **font)

def extract_data(filename):
    return np.fromregex(filename, r"\s*([^\s]+)\s+([^\s]+)", [('i', int), ('r2', float)])

def compute_D(r2):
    x = (np.arange(len(r2)) * dt)[:, None]
    y = r2 / 6

    # remove first skip instants
    x = x[skip:]
    y = y[skip:]

    model = LinearRegression(fit_intercept=False).fit(x, y)
    return model.coef_[0]

files = [
    # (filename, x)
    ("fort.8", 1)
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