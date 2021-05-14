import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


data = []
with open('100x100.txt', 'r') as file:
	for line in file:
		data.append([float(x) for x in line.split()])
	#print(len(data))
	#print(data)	 
	t = 0
	x = 0
	table = np.array(data)
	#np.resize(table, (100, 100))
	#table.append([])
	print(table)
	print(table.size)

	hf = plt.figure()
	ha = hf.add_subplot(111, projection='3d')
	x = np.linspace(0, 1, 101)
	y = np.linspace(0, 1, 101)	
	X, Y = np.meshgrid(x, y)
	ha.plot_surface(X, Y, table)
	plt.show()