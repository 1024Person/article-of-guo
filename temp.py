import numpy as np
from scipy.special import jv
import matplotlib.pyplot as plt

# v,z = np.indices([50,100])
# ax = plt.subplot(projection='3d')
# j=jv(v,z)
# ax.plot_surface(v,z,j)
# plt.show()
v = 1
z = -0.4
j = jv(v,z)
print(j)