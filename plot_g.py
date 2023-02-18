import matplotlib.pyplot as plt
import numpy as np

def extract_data(filename):
    return np.fromregex(filename, r"\s*([^\s]+)\s+([^\s]+)", [('rad', float), ('g', float)])

files = [
    # (filename, label)
]

plt.figure(figsize=(20,5))
for filename, label in files:
    data = extract_data(filename)
    r = data['rad']
    g = data['g']

    plt.plot(r, g, label=label)

plt.grid()
plt.title('g(r)')

if len(files) > 1:
    plt.legend()

plt.show()
