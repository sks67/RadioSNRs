import numpy as np

numbers = np.random.random(100)
chosen_bins1 = np.array([0.2,0.4,0.6,0.8])
chosen_bins2 = np.array([0.0,0.2,0.4,0.6,0.8,1.0])
n,bins = np.histogram(numbers, bins=chosen_bins1)
n2,bins2 = np.histogram(numbers, bins=chosen_bins2)
print n, np.sum(n)
print bins
print n2,np.sum(n2)
print bins2
binwidths = bins[1:]-bins[:-1]
print binwidths
