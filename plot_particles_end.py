import matplotlib.pyplot as plt
import numpy as np

with open('fort.1') as f:
    for line in f:
        pass
    last_line = line

data = last_line.split()
# remove iteration/time
data = data[2:]
data = data[:len(data) // 2]

data = np.array(tuple(map(float, data))).reshape(3, -1, order='F')
x = data[0]
y = data[1]
z = data[2]

ax = plt.figure(figsize=(20,5)).add_subplot(projection='3d')

ax.scatter(x,y,z)
ax.grid()

plt.show()
