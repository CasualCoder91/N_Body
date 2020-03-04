import matplotlib.pyplot as plt
import numpy as np
import imf

M = np.random.random(size=100000)
N = 0.08 + (0.2 * M**1.55 + 0.05 * M**0.6 ) / (1-M)**0.58

plt.hist(N, bins=np.logspace(-2, 2, 20), alpha=0.5, color="b")

masses = imf.make_cluster(10000, massfunc="kroupa")
plt.hist(masses, bins=np.logspace(-2, 2, 20), alpha=0.5, color="r")

plt.loglog ()
plt.show()