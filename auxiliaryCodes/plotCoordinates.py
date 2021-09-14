import matplotlib.pyplot as plt
import numpy as np

coordinatesLayer0 = np.loadtxt("../systemInfo/geometry/coordinatesLayer0")
coordinatesLayer1 = np.loadtxt("../systemInfo/geometry/coordinatesLayer1")

plt.scatter(coordinatesLayer0[:,0], coordinatesLayer0[:,1], c='red',  label='layer0')
plt.scatter(coordinatesLayer1[:,0], coordinatesLayer1[:,1], c='blue', label='layer1')
plt.legend()
plt.show()
