import numpy as np
import matplotlib.pyplot as plt

Presion=[]
Temperatura=[]

def my_range(start, end, step):
    while start <= end:
        yield start
        start += step

def imprimeDesv(i):
	filename = "presion{:{width}.{prec}f}".format(i, width=3, prec=2)
	filename = filename.replace('.',',')+".dat"
	file = open(filename,"r")
	velocidades = file.readlines()
	velocidades = list(map(float, velocidades))
	Presion.append(velocidades.pop())
	stde = np.std(velocidades)
	file.close()
	Temperatura.append(stde*stde)
	

for i in my_range(1.00,15.00,0.25):
	imprimeDesv(i)

for i in my_range(15.05,19.55,0.5):
	imprimeDesv(i)	
imprimeDesv(20)
imprimeDesv(30)

plt.plot(Temperatura,Presion,'.r')
plt.xlabel('Temperatura')
plt.ylabel('Presion')
plt.show()
