import matplotlib.pyplot as plt
import numpy as np
array = 10**(np.random.normal(loc=51.0,scale=0.28,size=1000000)-51.0)
bns = np.logspace(np.log10(array.min()),np.log10(array.max()),100)
plt.hist(array,bins=bns,histtype='step',color='k')
#plt.xlim(8,10.0)
#lt.ylim(0,100)
fig = plt.gcf()
fig.set_size_inches(6,6)
plt.xscale('log')
plt.show()
