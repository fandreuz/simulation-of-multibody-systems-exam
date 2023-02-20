import numpy as np

def extract(filename):
    with open(filename) as f:
        for line in f:
            pass
        last_line = line

    data = last_line.split()
    # remove iteration/time
    data = data[2:]
    data = data[:len(data) // 2]

    return np.array(tuple(map(float, data))).reshape(3, -1, order='F')