import numpy as np
import sys

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

def extract_to(filename, dest):
    data = extract(filename)
    np.savetxt(dest, data.T, fmt="%12.8f")

if __name__ == "__main__":
    extract_to(sys.argv[1], sys.argv[2])