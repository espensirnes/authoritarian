import numpy as np
from matplotlib import pyplot as plt



VOLATILITY = 0.01

#Defining the q function:
def q(pi, q0):
	ret = np.exp(-pi)+np.random.normal(scale=VOLATILITY)
	return ret



#This function runs the q function iteratively:
def path(its):

	p = []
	rho = 1
	q0 = 0.5
	pi0 = 0.9
	for i in range(its):
		q0 = q(pi0,q0)
		p.append(q0)

	return p

time = np.arange(100)

#Creating 50 different paths
for i in range(50):
	plt.plot(time,path(100))
plt.show()

a=0










