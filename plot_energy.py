import matplotlib.pyplot as plt
import numpy as np

data = np.fromregex('fort.2', r"\s*([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)", [('iteration', int), ('kinetic', float), ('potential', float), ('total', float), ('temperature', float)])
kinetic = data['kinetic']
potential = data['potential']
total = data['total']
temperature = data['temperature']

plt.figure(figsize=(20,5))

plt.subplot(1,2,1)
plt.plot(kinetic, label='Kinetic')
plt.plot(potential, label='Potential')
plt.plot(total, label='Total')
plt.grid()
plt.title('Energy')
plt.legend()

plt.subplot(1,2,2)
plt.plot(temperature)
plt.grid()
plt.title('Temperature')

plt.show()
