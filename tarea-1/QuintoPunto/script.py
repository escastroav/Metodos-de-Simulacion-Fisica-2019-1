import numpy as np
import matplotlib.pyplot as plt

Presion=[]
Temperatura=[]

def imprimeDesv(i):
	file = open("presion{}.dat".format(i),"r")
	velocidades = file.readlines()
	velocidades = list(map(float, velocidades))
	Presion.append(velocidades.pop())
	stde = np.std(velocidades)
	file.close()
	Temperatura.append(stde*stde)

for i in range(1.00,10.25,0.25):
	print(i)
	# imprimeDesv(i)
# imprimeDesv(15)
# imprimeDesv(20)
# imprimeDesv(30)

# plt.plot(Temperatura,Presion)
# plt.show()
