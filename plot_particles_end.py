import matplotlib.pyplot as plt
from extract_positions import extract

data = extract("fort.1")
x = data[0]
y = data[1]
z = data[2]

ax = plt.figure(figsize=(20,5)).add_subplot(projection='3d')

ax.scatter(x,y,z)
ax.grid()

plt.show()
